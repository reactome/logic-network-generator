import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from src.argument_parser import logger
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.logic_network_generator import (
    create_pathway_logic_network,
    export_cofactors,
    export_drugs,
    export_containment,
    export_entity_reaction_proxy_mapping,
    export_node_reaction_context,
    export_node_resolution,
    export_nodes,
    export_uuid_to_reactome_mapping,
)
from src.neo4j_connector import get_reaction_connections, get_reactome_release
from src.reaction_generator import get_decomposed_uid_mapping, prime_entity_caches

# Env vars that change what a generated network contains. Anything listed here
# is part of the cache fingerprint, so flipping one invalidates stale caches
# instead of silently reusing a network built under the old setting.
_FINGERPRINTED_ENV = (
    "LNG_MAX_VARIANTS",
    "LNG_COMPLEX_AS_NODE",
    "LNG_SET_EXPAND",
    "LNG_CATALYST_BUNDLE",
    "LNG_DIAGRAM_CONNECTIVITY",
    "LNG_DIAGRAM_BRIDGE",
    "LNG_DIAGRAM_DIR",
    "LNG_HANDOFF_EDGES",
    "LNG_HANDOFF_HUB_MAX",
    "LNG_SET_MEMBERS_OR",
    "LNG_DIAGRAM_SET_MEMBER",
    # Of these three, only LNG_EMIT_ONE_SIDED is a genuine cache gap: it is
    # read in reaction_generator.decompose_by_reactions, so it changes
    # best_matches.csv / decomposed_uid_mapping.csv, which ARE cached.
    # LNG_BOUNDARY_EXPANSION and LNG_COMPOSITION_EDGES are read inside
    # create_pathway_logic_network, which is re-run from scratch every time, so
    # listing them buys provenance in fingerprint.json at the cost of an
    # unnecessary re-fetch when they are flipped. Kept deliberately: knowing
    # which settings produced a catalog has been worth more than the re-fetch,
    # and src_sha256 already invalidates on any source change anyway.
    "LNG_BOUNDARY_EXPANSION",
    "LNG_BOUNDARY_HIERARCHY",
    "LNG_COMPOSITION_EDGES",
    "LNG_EMIT_ONE_SIDED",
    # Determinism controls: node ids are uuid4 and several selections iterate
    # sets, so hash seeding changes emitted content (~5.8% of TP53 edges per
    # bin/create-pathways.py). A cache built unseeded is not comparable to a
    # seeded one.
    "PYTHONHASHSEED",
    "LNG_ALLOW_NONDETERMINISM",
)

_FINGERPRINT_FILE = "fingerprint.json"


def _cache_fingerprint() -> Dict[str, Any]:
    """Identity of everything that determines a cached network's contents.

    The per-pathway `cache/` dir is reused on file existence alone, with no
    record of what produced it. That makes a benchmark delta unattributable: a
    rerun with changed code or a changed LNG_* setting silently replays the old
    network (verified: the same env produced 139 edges from cache vs 58 fresh).
    Hashing the source tree catches uncommitted edits, which a git rev would not.
    """
    digest = hashlib.sha256()
    sources = sorted((Path(__file__).resolve().parent).glob("*.py"))
    # The CLI entrypoint also shapes output (it sets PYTHONHASHSEED and parses
    # the pathway list), so include it when present.
    entrypoint = Path(__file__).resolve().parents[1] / "bin" / "create-pathways.py"
    if entrypoint.exists():
        sources.append(entrypoint)
    for path in sources:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return {
        "src_sha256": digest.hexdigest(),
        "env": {key: os.environ.get(key, "") for key in _FINGERPRINTED_ENV},
        # Two different databases can both report release 97, so the release
        # number alone does not identify the source graph. Store a HASH of the
        # connection URL, never the URL itself — it may embed credentials, and
        # this file is written into the output tree.
        "db_identity": hashlib.sha256(
            os.environ.get("NEO4J_URL", "bolt://localhost:7687").encode()
        ).hexdigest()[:16],
        "reactome_release": _reactome_release_cached(),
    }


_RELEASE_CACHE: Dict[str, Any] = {}


def _reactome_release_cached() -> Any:
    """`get_reactome_release()` memoized for the process (a ~33 ms DB query)."""
    if "value" not in _RELEASE_CACHE:
        _RELEASE_CACHE["value"] = get_reactome_release()
    return _RELEASE_CACHE["value"]


def _cache_is_reusable(cache_dir: Path, current: Dict[str, Any]) -> bool:
    """True when `cache_dir` was produced by the current code/env/release.

    A cache with no fingerprint is *adopted* (stamped and reused) so existing
    output trees keep working without a forced full regeneration; protection
    starts from the next run. A cache whose fingerprint differs is rejected and
    regenerated, naming what changed.
    """
    path = cache_dir / _FINGERPRINT_FILE
    if not path.exists():
        if not any(cache_dir.glob("*.csv")):
            # Nothing cached yet — this is a fresh generation, not an adoption.
            # Warning here would report "reusing" an empty directory, and
            # stamping "adopted" would certify the cache this run is about to
            # produce as being of unverified provenance.
            return False

        # Adopt rather than force a full regeneration of existing output trees,
        # but do NOT stamp the current source hash onto a cache of unknown
        # provenance — that would certify it as "built by this code" when we
        # have no idea what built it. Record the adoption instead, so a later
        # run can still tell the provenance was never verified.
        logger.warning(
            "Cache %s has no fingerprint: reusing it, but its provenance is "
            "UNVERIFIED (it may have been built by different code). Delete the "
            "cache dir to force a clean, attributable regeneration.",
            cache_dir,
        )
        _write_cache_fingerprint(cache_dir, current, provenance="adopted")
        return True

    try:
        stored = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        logger.warning("Unreadable cache fingerprint %s (%s); regenerating.", path, exc)
        return False

    if stored == current:
        return True

    changed = []
    adopted = stored.get("provenance") == "adopted"
    if adopted:
        logger.warning(
            "Cache %s was adopted without provenance; source changes since then "
            "cannot be detected.",
            cache_dir,
        )
    elif stored.get("src_sha256") != current["src_sha256"]:
        changed.append("generator source")
    if stored.get("db_identity") != current["db_identity"]:
        changed.append("Neo4j connection target")
    # An UNKNOWN release is not evidence of a change. get_reactome_release()
    # returns None when the DBInfo lookup fails, so comparing None against a
    # stored 97 would force a spurious full regeneration on a transient blip.
    # Only treat the release as differing when both sides actually know it.
    stored_release = stored.get("reactome_release")
    current_release = current["reactome_release"]
    if (
        stored_release is not None
        and current_release is not None
        and stored_release != current_release
    ):
        changed.append(f"Reactome release ({stored_release} -> {current_release})")
    for key in _FINGERPRINTED_ENV:
        old = (stored.get("env") or {}).get(key, "")
        new = current["env"][key]
        if old != new:
            changed.append(f"{key} ({old!r} -> {new!r})")

    if not changed:
        # The only difference was an unknown release; keep the cache and leave
        # the existing stamp in place rather than discarding work.
        logger.info(
            "Cache fingerprint in %s differs only in an unknown Reactome release; "
            "reusing the cache.",
            cache_dir,
        )
        return True
    logger.warning(
        "Cache in %s was built under different inputs (%s); regenerating so the "
        "result reflects the current configuration.",
        cache_dir,
        "; ".join(changed) or "unknown difference",
    )
    return False


def _write_cache_fingerprint(
    cache_dir: Path, fingerprint: Dict[str, Any], provenance: str = "generated"
) -> None:
    record = dict(fingerprint)
    record["provenance"] = provenance
    if provenance == "adopted":
        # We did not build these files; claiming a source hash would be a lie.
        record["src_sha256"] = None
    try:
        (cache_dir / _FINGERPRINT_FILE).write_text(json.dumps(record, indent=2, sort_keys=True))
    except OSError as exc:
        logger.warning("Could not write cache fingerprint: %s", exc)



def generate_pathway_file(
    pathway_id: str,
    pathway_name: str,
    output_dir: str = "output",
) -> None:
    """Generate pathway logic network file with caching.

    Args:
        pathway_id: Reactome pathway database ID
        pathway_name: Human-readable pathway name
        output_dir: Base output directory (default: "output")

    Raises:
        ConnectionError: If Neo4j database is not accessible
        ValueError: If pathway data is invalid or pathway not found
        IOError: If cache files cannot be written

    Output files are organized as:
        {output_dir}/{pathway_name}_{pathway_id}/
            logic_network.csv           - Main logic network (what users need)
            stid_to_uuid_mapping.csv    - Stable ID to UUID mapping (what users need)
            cache/                      - Intermediate files
    """
    logger.info(f"Generating logic network for pathway {pathway_id}: {pathway_name}")

    # Create pathway-specific output directory
    base_output_dir = Path(output_dir)
    base_output_dir.mkdir(exist_ok=True)

    # Folder name is the stable ID and nothing else.
    #
    # It used to be "{sanitized_pathway_name}_{pathway_id}", which put a NAME in
    # an identifier. Names are data: they gain commas, lose trailing
    # underscores, and get recurated. On Release97 three differ between the
    # curator files and these directories -- "Interleukin-3,_Interleukin-5...",
    # "Signaling_by_..._IGF1R_", "Mitotic_G1-G1_S_phases" -- and matching on
    # name silently misfiled 438 of 1,484 cases in an analysis whose cases all
    # came from this catalog.
    #
    # The pathway name is still recoverable: it is in the Reactome database
    # under this id, and consumers that want it should look it up rather than
    # parse it out of a path. Downstream lookups already key on the id suffix.
    # Guard before the id becomes a path. Path("out") / "" is "out", so an empty
    # id would make the pathway directory the CATALOG ROOT and write
    # logic_network.csv, cache/ and the rest straight into it, colliding with
    # every other pathway. The old "{name}_{id}" spelling hid this because the
    # name kept the directory distinct; naming by id alone makes it reachable.
    folder_name = str(pathway_id).strip()
    if not folder_name or folder_name in {".", "..", "None"} or "/" in folder_name:
        raise ValueError(
            f"Refusing to generate with an unusable pathway id {pathway_id!r}: "
            "the directory is named by the id, so this would write into the "
            "catalog root."
        )
    pathway_output_dir = base_output_dir / folder_name
    pathway_output_dir.mkdir(exist_ok=True)

    # Create cache subdirectory for intermediate files
    cache_dir = pathway_output_dir / "cache"
    cache_dir.mkdir(exist_ok=True)

    # Define filenames for caching (in cache subdirectory)
    reaction_connections_file = cache_dir / "reaction_connections.csv"
    decomposed_uid_mapping_file = cache_dir / "decomposed_uid_mapping.csv"
    best_matches_file = cache_dir / "best_matches.csv"

    try:
        # Decide ONCE per pathway whether this cache dir may be reused. Without
        # this the reuse decision is a bare exists() check with no record of what
        # produced the files, so a rerun after a code or LNG_* change silently
        # replays the old network and a benchmark delta cannot be attributed.
        fingerprint = _cache_fingerprint()
        cache_usable = _cache_is_reusable(cache_dir, fingerprint)

        # Load or fetch reaction connections
        if cache_usable and os.path.exists(reaction_connections_file):
            logger.info(f"Loading cached reaction connections from {reaction_connections_file}")
            reaction_connections = pd.read_csv(reaction_connections_file, dtype=str)
        else:
            logger.info(f"Fetching reaction connections from Neo4j for pathway {pathway_id}")
            reaction_connections = get_reaction_connections(pathway_id)
            try:
                reaction_connections.to_csv(reaction_connections_file, index=False)
                logger.info(f"Cached reaction connections to {reaction_connections_file}")
            except IOError as e:
                logger.warning(f"Could not cache reaction connections: {e}")
                # Continue without caching

        # Load or generate decomposition and best matches
        if cache_usable and os.path.exists(decomposed_uid_mapping_file) and os.path.exists(best_matches_file):
            logger.info(f"Loading cached decomposition from {decomposed_uid_mapping_file}")
            decomposed_uid_mapping = pd.read_csv(
                decomposed_uid_mapping_file,
                dtype=decomposed_uid_mapping_column_types,  # type: ignore[arg-type]
            )
            best_matches = pd.read_csv(best_matches_file)
            # Reusing the cached decomposition skips get_decomposed_uid_mapping,
            # which is the only place the entity caches are cleared and primed.
            # Prime them here or this pathway inherits the previous pathway's
            # prefetch state and every one of its complexes resolves as atomic.
            prime_entity_caches(reaction_connections)
        else:
            logger.info("Decomposing complexes and entity sets...")
            [decomposed_uid_mapping, best_matches_list] = get_decomposed_uid_mapping(
                pathway_id, reaction_connections
            )
            best_matches = pd.DataFrame(
                best_matches_list,
                columns=["incomming", "outgoing", "reactome_id"],
            )

            try:
                decomposed_uid_mapping.to_csv(decomposed_uid_mapping_file, index=False)
                best_matches.to_csv(best_matches_file, index=False)
                logger.info(f"Cached decomposition to {decomposed_uid_mapping_file}")
                # Stamp the cache we just PRODUCED with real provenance. Without
                # a write here the mechanism only ever adopts: generate, change
                # the code, rerun, and the next run would adopt the stale cache
                # and stamp the NEW hash onto it.
                _write_cache_fingerprint(cache_dir, fingerprint, provenance="generated")
            except IOError as e:
                logger.warning(f"Could not cache decomposition results: {e}")
                # Continue without caching

        # Add curator-drawn diagram connectivity (product->substrate pairs the
        # diagram links but precedingEvent may omit). See reactome/logic-network-generator#39.
        #
        # DEFAULT = MERGE (LNG_DIAGRAM_BRIDGE=0): a diagram-drawn connection is
        # made the same way an annotated one is — by merging the shared product
        # into a single node, so A -> product -> B. One drawn connection becomes
        # one connection.
        #
        # This is needed and it is not optional. Of 650 human pathways that have
        # their own diagram with reaction glyphs, 89 (13.7%) have under half
        # their reactions carrying a precedingEvent. Without diagram
        # connectivity those come out as disconnected node piles.
        #
        # It used to default to ADDITIVE bridges, which connected the producer's
        # output-copy to the consumer's input-copy as new edges. That is a
        # cartesian product over variant instances: |variants of A| x |variants
        # of B| x |shared entities|. Measured, it averaged ~35 edges per drawn
        # connection and 577 per connection that actually needed bridging, which
        # is how diagram_bridge grew to 14.3% of every edge in the catalog. On
        # the wide curator set, held out from tuning, removing those edges was
        # worth +107 cases; merging instead is +1, i.e. free.
        #
        # Deliberately NOT gated on whether the consumer reaction already has a
        # precedingEvent. Curators annotate incrementally, so a reaction with
        # one annotated predecessor may still be missing others — presence of
        # annotation does not mean completeness. And absence does not mean a
        # gap: of six pathways with 0% coverage, four have no chainable reaction
        # pairs at all, so their emptiness is correct. Neither signal is
        # reliable, so we do not guess: every drawn connection is made, once.
        #
        # Set LNG_DIAGRAM_BRIDGE=1 for the legacy additive bridges;
        # set LNG_DIAGRAM_CONNECTIVITY=0 to disable diagram connectivity entirely.
        from src.diagram_connectivity import (
            augment_reaction_connections,
            diagram_set_member_pairs,
            diagram_shared_product_pairs,
        )
        diagram_bridge_pairs = None
        if os.environ.get("LNG_DIAGRAM_BRIDGE", "0") == "1":
            connectivity = reaction_connections
            diagram_bridge_pairs = diagram_shared_product_pairs(pathway_id)
        else:
            connectivity = augment_reaction_connections(pathway_id, reaction_connections)

        # Generate logic network
        logger.info("Creating pathway logic network...")
        # Realisation links the diagram draws between a specific complex and
        # the generic one it instantiates. Has its own flag AND honours
        # LNG_DIAGRAM_CONNECTIVITY=0, so the documented "disable diagram
        # connectivity entirely" kill switch above stays true.
        # DEFAULT OFF. The edge is correct in principle and inert in practice:
        # the target complex already carries `and` assembly edges to every
        # constituent protein, and DeltaSignal combines an AND cluster with an
        # OR cluster as max(and, or), so a DECREASE through this edge is
        # discarded — measured on TGF-beta, member DOWN leaves the generic at
        # 1.0 while member UP reaches 80.0. One-directional propagation is
        # worse than none, because it looks like it works.
        #
        # It becomes load-bearing under DS_OR_COMBINE=gate, which is itself
        # default-off and measured inert precisely because nothing in the
        # catalog produced mixed and/or clusters. This produces them. The two
        # therefore have to be enabled and measured TOGETHER, and neither alone.
        set_member_pairs = None
        if (os.environ.get("LNG_DIAGRAM_SET_MEMBER", "0") == "1"
                and os.environ.get("LNG_DIAGRAM_CONNECTIVITY", "1") != "0"):
            try:
                set_member_pairs = diagram_set_member_pairs(pathway_id)
            except Exception:
                logger.warning("Could not read diagram set-member links",
                               exc_info=True)

        result = create_pathway_logic_network(
            decomposed_uid_mapping, connectivity, best_matches,
            diagram_bridge_pairs=diagram_bridge_pairs,
            diagram_set_member_pairs=set_member_pairs,
        )

        # Save logic network (main output file users need)
        output_file = pathway_output_dir / "logic_network.csv"
        try:
            result.logic_network.to_csv(output_file, index=False)
            logger.info(f"Successfully generated logic network: {output_file}")
            logger.info(f"Network contains {len(result.logic_network)} edges")
        except IOError as e:
            logger.error(f"Failed to write output file {output_file}: {e}")
            raise

        # Export UUID to Reactome stable ID mapping (main mapping file users need)
        uuid_to_reactome_file = pathway_output_dir / "stid_to_uuid_mapping.csv"
        try:
            export_uuid_to_reactome_mapping(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                result.catalyst_regulator_map,
                str(uuid_to_reactome_file)
            )
            logger.info(f"Successfully exported stable ID to UUID mapping: {uuid_to_reactome_file}")
        except IOError as e:
            logger.error(f"Failed to write stable ID to UUID mapping file {uuid_to_reactome_file}: {e}")
            # Don't raise - this is supplementary

        # Export entity→reaction proxy mapping. Curated species (often Complexes)
        # that were expanded into virtual variants lose their own stId from the
        # UUID mapping; this file points each such species at the UUIDs of the
        # reactions that produce (or, failing that, consume) it, so consumers can
        # read reaction flux as a proxy for the species' state.
        proxy_mapping_file = pathway_output_dir / "entity_reaction_proxy_mapping.csv"
        try:
            export_entity_reaction_proxy_mapping(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                pathway_id,
                str(proxy_mapping_file),
            )
            logger.info(f"Successfully exported entity-reaction proxy mapping: {proxy_mapping_file}")
        except Exception as e:
            logger.error(f"Failed to write entity-reaction proxy mapping file {proxy_mapping_file}: {e}")
            # Don't raise - this is supplementary

        # Schema-backed provenance files (schema/logic_network.linkml.yaml):
        # nodes.csv (node_kind, diagram_entity_id, member_leaves, set provenance)
        # and node_reaction_context.csv (node<->reaction location layer).
        try:
            export_nodes(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                str(pathway_output_dir / "nodes.csv"),
            )
            export_node_reaction_context(
                result.entity_uuid_registry,
                result.reaction_id_map,
                result.catalyst_regulator_map,
                str(pathway_output_dir / "node_reaction_context.csv"),
                logic_network=result.logic_network,
            )
            export_node_resolution(
                pathway_id,
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                str(pathway_output_dir / "node_resolution.csv"),
                str(pathway_output_dir / "node_exclusions.csv"),
            )
            # Ships WITH the networks so an artifact bundle pulled from S3
            # carries its own answer to "which of these nodes is ATP".
            export_cofactors(
                result.logic_network,
                result.uuid_mapping,
                str(pathway_output_dir / "cofactors.csv"),
            )

            # Drug-derived entities (deltasignal specs/032), for a consumer that
            # models a cell without the drug.
            export_drugs(
                result.logic_network,
                result.uuid_mapping,
                str(pathway_output_dir / "drugs.csv"),
            )

            # What each node CONTAINS, so a consumer can select nodes by
            # containment instead of us inventing assembly/dissociation edges.
            export_containment(
                result.uuid_mapping,
                str(pathway_output_dir / "containment.csv"),
            )
        except Exception as e:
            logger.error(f"Failed to write node provenance files: {e}", exc_info=True)
            # Don't raise - supplementary

        logger.info(f"Output directory: {pathway_output_dir}")

    except (ConnectionError, ValueError) as e:
        logger.error(f"Failed to generate pathway {pathway_id}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error generating pathway {pathway_id}", exc_info=True)
        raise RuntimeError(f"Pathway generation failed: {str(e)}") from e
