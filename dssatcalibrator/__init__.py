"""dssatcalibrator — Monte-Carlo & Bayesian calibration of DSSAT-CSM.

Concept implementation. See ../CONCEPT.md for the architecture and decisions.

Public surface is intentionally small and grows as modules land:
    dssat_io      — parse DSSAT outputs (PlantGro.OUT, Evaluate.OUT, Summary.OUT)
    observations  — read observed data (FileA/FileT, long-format CSV)
"""

__version__ = "0.0.1"
