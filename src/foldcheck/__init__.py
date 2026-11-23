"""Score AlphaFold models for geometry problems and estimated accuracy."""

__version__ = "0.1.0"

from foldcheck.scoring import score_file, score_structure

__all__ = ["__version__", "score_file", "score_structure"]
