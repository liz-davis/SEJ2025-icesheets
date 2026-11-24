"""Scientific constants preserved from the SEJ2025 production chain.

No paths are created on import. All source files are supplied at runtime.
"""

from pathlib import Path
import numpy as np

YEARS = np.arange(2000, 2301)
SCENARIOS = ("LOW", "INTERMEDIATE", "HIGH")
SCENARIO_CODES = {"LOW": "2", "INTERMEDIATE": "35", "HIGH": "5"}
REFERENCE_YEARS = np.array([2000, 2020, 2050, 2100, 2200, 2300])
REFERENCE_WARMING = {
    "LOW": np.array([0.7, 1.2, 1.7, 2.0, 2.0, 2.0]),
    "INTERMEDIATE": np.array([0.7, 1.2, 2.0, 3.5, 3.5, 3.5]),
    "HIGH": np.array([0.7, 1.2, 2.5, 5.0, 5.0, 5.0]),
}
CHECKPOINTS = (2100, 2200)
SELECTION_SEED = 8172026
CMIP6_SAMPLE_SEED = 83126
CMIP7_CRN_SEED = 100626
TENSION_P = 0.03
COMPONENTS = ("E", "W", "G")
DEFAULT_COUPLED_FILENAME = "SEJ25_emulandice2050_recoupled_by_warming.csv"


def repository_root() -> Path:
    """Repo root for editable src/ installations; no dependence on cwd."""
    return Path(__file__).resolve().parents[2]
