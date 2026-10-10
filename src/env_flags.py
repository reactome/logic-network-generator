"""One reader for the generator's LNG_* switches.

Before this (code review 2026-10-02) each flag was compared to a literal at its
call site, in two opposite ways: `!= "0"` flags turned ON for any other value
(LNG_HANDOFF_EDGES=false enabled the measured-harmful handoff edges) and
`== "1"` flags turned OFF (LNG_BOUNDARY_HIERARCHY=true disabled the adopted
hierarchy). A misspelt name was simply never read, while catalog.sh recorded it
in BUILD.json as if applied. Now a switch takes exactly "0" or "1", and an
LNG_* name the generator does not read is an error at startup.
"""
import os

# Boolean switches and their defaults.
BOOL_FLAGS = {
    "LNG_BOUNDARY_EXPANSION": True,
    "LNG_BOUNDARY_HIERARCHY": True,
    "LNG_BIND_STOICH": True,     # default since deltasignal specs/048 (homodimer binding; 0 predictions moved)
    "LNG_CAP_POOLS": True,      # default since deltasignal specs/045 (RAF experimental +14, held-out 0)
    "LNG_CATALYST_BUNDLE": False,
    "LNG_COMPLEX_AS_NODE": True,
    "LNG_COMPOSITION_EDGES": False,
    "LNG_DEBUG_TRACEBACKS": False,
    "LNG_DIAGRAM_BRIDGE": False,
    "LNG_DIAGRAM_CONNECTIVITY": True,
    "LNG_DIAGRAM_SET_MEMBER": False,
    "LNG_EMIT_ONE_SIDED": True,
    "LNG_HANDOFF_EDGES": False,
    "LNG_SET_EXPAND": True,
    "LNG_SET_MEMBERS_OR": False,
    "LNG_SET_POOL": True,
    "LNG_VARIANT_NODES": True,   # default since deltasignal specs/046 (vn7: held-out +130, experimental -10 ns)
}

# Read elsewhere with their own parsing and validation.
OTHER_FLAGS = frozenset({
    "LNG_ALLOW_NONDETERMINISM",   # bin/create-pathways.py, before imports
    "LNG_DIAGRAM_DIR",            # a path
    "LNG_HANDOFF_HUB_MAX",        # _int_env
    "LNG_MAX_VARIANTS",           # _int_env
    "LNG_POOL_ACTIVE_VIA",        # pool_active_via()
    "LNG_PYTHON",                 # catalog.sh: the interpreter, not a switch
})


def env_flag(name: str) -> bool:
    """The value of boolean switch `name`: unset or empty is its default,
    "1" is on, "0" is off, anything else is an error."""
    default = BOOL_FLAGS[name]          # KeyError: an unregistered switch
    raw = os.environ.get(name, "").strip()
    if raw == "":
        return default
    if raw not in ("0", "1"):
        raise ValueError(f"{name} must be 0 or 1, got {raw!r}. Unset it for the default "
                         f"({int(default)}).")
    return raw == "1"


def validate_env(removed=()) -> None:
    """Every LNG_* in the environment is a known name with a valid value."""
    known = set(BOOL_FLAGS) | OTHER_FLAGS | set(removed)
    unknown = sorted(n for n in os.environ if n.startswith("LNG_") and n not in known)
    if unknown:
        raise ValueError(f"Unknown LNG_* variable(s) {unknown}: the generator does not read "
                         "them, so they would have no effect. Check the spelling.")
    for name in BOOL_FLAGS:
        env_flag(name)
