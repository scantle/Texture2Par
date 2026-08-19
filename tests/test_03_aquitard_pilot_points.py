"""Regression coverage for independent aquifer/aquitard pilot-point zones."""

from __future__ import annotations

import csv


DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,1.0,1.0,100.0,0.0,2.0,0.6,0.4
W1,1,2,1.0,1.0,100.0,2.0,22.0,0.6,0.4
W2,2,1,5.0,5.0,100.0,0.0,2.0,0.6,0.4
W2,2,2,5.0,5.0,100.0,2.0,22.0,0.6,0.4
W3,3,1,9.0,1.0,100.0,0.0,2.0,0.6,0.4
W3,3,2,9.0,1.0,100.0,2.0,22.0,0.6,0.4
"""

PREPROCESSOR = """\
Title 1
Title 2
Title 3
preprocessor.bin
elements.dat
nodes.dat
stratigraphy.dat
"""

SIMULATION = """\
Title 1
Title 2
Title 3
preprocessor.bin
groundwater.dat
"""

NODES = """\
6
1.0
1 0.0 0.0
2 5.0 0.0
3 10.0 0.0
4 0.0 10.0
5 5.0 10.0
6 10.0 10.0
"""

ELEMENTS = """\
4
1
1 TestRegion
1 1 2 4 0
2 2 5 4 0
3 2 3 5 0
4 3 6 5 0
"""

STRATIGRAPHY = """\
1
1.0
1 100.0 2.0 20.0
2 100.0 2.0 20.0
3 100.0 2.0 20.0
4 100.0 2.0 20.0
5 100.0 2.0 20.0
6 100.0 2.0 20.0
"""

PP_ZONES = """\
2 1
# Aquifer array
CONSTANT 1
# Aquitard array
CONSTANT 1
"""

COMMON_FILES = {
    "dataset.csv": DATASET,
    "preprocessor.in": PREPROCESSOR,
    "simulation.in": SIMULATION,
    "nodes.dat": NODES,
    "elements.dat": ELEMENTS,
    "stratigraphy.dat": STRATIGRAPHY,
}


def _input(*, explicit_zones: bool) -> str:
    zone_file = "PP_ZONE_FILE ppzones.dat" if explicit_zones else ""
    aquifer_zone = 1 if explicit_zones else 7
    aquitard_zone = 1 if explicit_zones else 2
    return f"""\
BEGIN OPTIONS
  MAX_VSTRUCT 1
  INFER_LAST_CLASS
  WRITE_NODE_FILES
  WRITE_NONE
END OPTIONS

BEGIN FLOW_MODEL
  TYPE IWFM
  SIM_FILE simulation.in
  PREPROC_FILE preprocessor.in
  {zone_file}
END FLOW_MODEL

BEGIN CLASSES
  Coarse
  Fine
END CLASSES

BEGIN DATASET
  FILE dataset.csv
END DATASET

BEGIN VARIOGRAMS
  CLASS Coarse
    1 Sph 0.0 1.0 100.0 100.0 0.0 1
  CLASS PilotPoints
    1 Sph 0.0 1.0 100.0 100.0 0.0 1
END VARIOGRAMS

BEGIN PP_LOCS
  1 1.0 1.0 {aquifer_zone}
  2 5.0 5.0 {aquifer_zone}
  3 9.0 1.0 {aquifer_zone}
END PP_LOCS

BEGIN PP_AQUITARD_LOCS
  1 1.0 1.0 {aquitard_zone}
  2 5.0 5.0 {aquitard_zone}
  3 9.0 1.0 {aquitard_zone}
END PP_AQUITARD_LOCS

BEGIN PP_PARAMETERS
  TYPE GLOBAL
  1 1.0 1.0 1.0
  2 1.0 1.0 1.0
  3 1.0 1.0 1.0

  TYPE AQUIFER
  1 Coarse 100.0 100.0 1e-5 0.15 1.0 0.0
  1 Fine   100.0 100.0 1e-5 0.15 1.0 0.0
  2 Coarse 100.0 100.0 1e-5 0.15 1.0 0.0
  2 Fine   100.0 100.0 1e-5 0.15 1.0 0.0
  3 Coarse 100.0 100.0 1e-5 0.15 1.0 0.0
  3 Fine   100.0 100.0 1e-5 0.15 1.0 0.0

  TYPE GLOBAL_AQUITARD
  1 1.0 1.0
  2 1.0 1.0
  3 1.0 1.0

  TYPE AQUITARD
  1 Coarse 0.01 0.01 1.0 0.0
  1 Fine   0.01 0.01 1.0 0.0
  2 Coarse 0.01 0.01 1.0 0.0
  2 Fine   0.01 0.01 1.0 0.0
  3 Coarse 0.01 0.01 1.0 0.0
  3 Fine   0.01 0.01 1.0 0.0
END PP_PARAMETERS
"""


def _column_values(path, column):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return [float(row[column]) for row in rows]


def test_default_pp_zones_match_explicit_constant_zones(run_t2p):
    """No zone file means zone 1 independently for aquifer and aquitard PPs."""
    no_zone = run_t2p(
        _input(explicit_zones=False),
        filename="no_zone.t2p",
        extra_files=COMMON_FILES,
    )
    assert no_zone.returncode == 0, no_zone.combined
    no_zone_aquifer = no_zone.output_file("no_zone_KV.csv")
    no_zone_aquitard = no_zone.output_file("no_zone_KV_aqtard.csv")
    assert no_zone_aquifer.exists()
    assert no_zone_aquitard.exists()
    no_zone_aquifer_text = no_zone_aquifer.read_text(encoding="utf-8")
    no_zone_aquitard_text = no_zone_aquitard.read_text(encoding="utf-8")
    no_zone_aquifer_kv = _column_values(no_zone_aquifer, "Layer1")
    no_zone_aquitard_kv = _column_values(no_zone_aquitard, "Aquitard1")

    with_zone = run_t2p(
        _input(explicit_zones=True),
        filename="with_zone.t2p",
        extra_files={**COMMON_FILES, "ppzones.dat": PP_ZONES},
    )
    assert with_zone.returncode == 0, with_zone.combined
    with_zone_aquifer = with_zone.output_file("with_zone_KV.csv")
    with_zone_aquitard = with_zone.output_file("with_zone_KV_aqtard.csv")
    assert with_zone_aquifer.exists()
    assert with_zone_aquitard.exists()

    assert with_zone_aquifer.read_text(encoding="utf-8") == no_zone_aquifer_text
    assert with_zone_aquitard.read_text(encoding="utf-8") == no_zone_aquitard_text
    assert _column_values(with_zone_aquifer, "Layer1") == no_zone_aquifer_kv
    assert _column_values(with_zone_aquitard, "Aquitard1") == no_zone_aquitard_kv
    assert all(abs(value - 100.0) < 1e-10 for value in no_zone_aquifer_kv)
    assert all(abs(value - 0.01) < 1e-10 for value in no_zone_aquitard_kv)
