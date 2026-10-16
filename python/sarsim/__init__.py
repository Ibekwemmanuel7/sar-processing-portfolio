"""Stripmap SAR point-target simulator, range-Doppler processor and verification tools.

Separate Python implementation of the MATLAB code in ../matlab, used for numerical comparison.
It shares the same parameters and model assumptions, so it is not an independent validation.
"""
from .core import (params, simulate_raw, rda_focus, measure_irf, irf_model, nesz_model,
                   slick_scene, detect_dark)

__all__ = ["params", "simulate_raw", "rda_focus", "measure_irf", "irf_model", "nesz_model",
           "slick_scene", "detect_dark"]
