"""CAPPED_IDS survives a cached rebuild and does not leak between pathways
(code review F6; LNG_CAP_POOLS default since deltasignal specs/045).

CAPPED_IDS is filled only while decomposing. A rebuild that reused the cached
decomposition skipped that, so the network builder saw no capped reaction and
emitted no cap pools; and the set was never cleared, so ids leaked from one
pathway to the next. Driven through generate_pathway_file."""
import pandas as pd
import pytest

import src.pathway_generator as pg
import src.reaction_generator as rg
from tests.test_export_failfast import mocked  # noqa: F401  (fixture)


@pytest.fixture
def seen(mocked):  # noqa: F811
    record = []
    mocked.setattr(pg, "_cache_fingerprint", lambda: {"fixed": 1})
    mocked.setattr(pg, "_write_cache_fingerprint", lambda d, fp, provenance=None:
                   (d / "fingerprint.json").write_text("{}"))
    mocked.setattr(pg, "_cache_is_reusable", lambda d, fp: (d / "fingerprint.json").exists())
    mocked.setattr(pg, "prime_entity_caches", lambda rc: None)
    real = pg.create_pathway_logic_network

    def spy(*a, **k):
        record.append(set(rg.CAPPED_IDS))
        return real(*a, **k)
    mocked.setattr(pg, "create_pathway_logic_network", spy)
    return mocked, record


def decompose_capping(ids):
    def f(pid, rc):
        rg.CAPPED_IDS.update(ids)
        return pd.DataFrame({"uid": ["u"]}), [("a", "b", "R-HSA-1")]
    return f


def test_a_cached_rebuild_restores_the_capped_ids(tmp_path, seen):
    mp, record = seen
    mp.setattr(pg, "get_decomposed_uid_mapping", decompose_capping({"R-HSA-CAP"}))
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    assert (tmp_path / "R-HSA-1" / "cache" / "capped_ids.txt").read_text() == "R-HSA-CAP\n"

    def must_not_decompose(pid, rc):
        raise AssertionError("the cache should have been reused")
    mp.setattr(pg, "get_decomposed_uid_mapping", must_not_decompose)
    rg.CAPPED_IDS.clear()
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    assert record == [{"R-HSA-CAP"}, {"R-HSA-CAP"}]


def test_capped_ids_do_not_leak_between_pathways(tmp_path, seen):
    mp, record = seen
    mp.setattr(pg, "get_decomposed_uid_mapping", decompose_capping({"R-HSA-CAP"}))
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    mp.setattr(pg, "get_decomposed_uid_mapping", decompose_capping(set()))
    pg.generate_pathway_file("R-HSA-2", "y", str(tmp_path))
    assert record == [{"R-HSA-CAP"}, set()]


def test_a_cache_without_the_capped_ids_file_is_not_reused(tmp_path, seen):
    mp, record = seen
    mp.setattr(pg, "get_decomposed_uid_mapping", decompose_capping({"R-HSA-CAP"}))
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    (tmp_path / "R-HSA-1" / "cache" / "capped_ids.txt").unlink()
    calls = []

    def again(pid, rc):
        calls.append(1)
        rg.CAPPED_IDS.update({"R-HSA-CAP"})
        return pd.DataFrame({"uid": ["u"]}), [("a", "b", "R-HSA-1")]
    mp.setattr(pg, "get_decomposed_uid_mapping", again)
    pg.generate_pathway_file("R-HSA-1", "x", str(tmp_path))
    assert calls == [1] and record[-1] == {"R-HSA-CAP"}
