"""Stripmap SAR point-target simulator, range-Doppler processor and verification tools.

Python port of the MATLAB code in ../matlab (same parameters and algorithms).
"""
from .core import (params, simulate_raw, rda_focus, measure_irf, irf_model, nesz_model,
                   slick_scene, detect_dark)

__all__ = ["params", "simulate_raw", "rda_focus", "measure_irf", "irf_model", "nesz_model",
           "slick_scene", "detect_dark"]
