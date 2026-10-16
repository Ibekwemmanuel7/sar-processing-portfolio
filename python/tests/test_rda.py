"""Regression tests: simulator + range-Doppler processor meet the requirement and match the model."""
import copy

import numpy as np
import pytest

from sarsim import params, simulate_raw, rda_focus, measure_irf, irf_model, nesz_model, slick_scene, detect_dark


@pytest.fixture(scope="module")
def single():
    p = params(); p.targets = np.array([[0, 0, 1.0]]); p.snr_in_db = 40
    raw, t, eta = simulate_raw(p)
    return p, raw, t, eta


def test_weighted_meets_requirement(single):
    p, raw, t, eta = single
    img, _ = rda_focus(raw, t, p)
    m = measure_irf(img, p, t, eta, 0, 0); mdl = irf_model(p); s = p.spec
    assert m.irw_rg_m <= s.irw_rg_m and m.irw_az_m <= s.irw_az_m
    assert max(m.pslr_rg_db, m.pslr_az_db) <= s.pslr_db
    assert max(m.islr_rg_db, m.islr_az_db) <= s.islr_db
    assert abs(m.irw_rg_m / mdl.irw_rg_m - 1) <= s.model_tol
    assert abs(m.irw_az_m / mdl.irw_az_m - 1) <= s.model_tol
    assert max(abs(m.loc_err_rg_m), abs(m.loc_err_az_m)) <= s.loc_err_m


def test_unweighted_sinc_sidelobes(single):
    p, raw, t, eta = single
    q = copy.copy(p); q.win_alpha = 1.0
    img, _ = rda_focus(raw, t, q)
    m = measure_irf(img, q, t, eta, 0, 0)
    assert abs(m.pslr_rg_db + 13.26) < 1.0


def test_nesz_requirement():
    p = params(); _, _, nesz = nesz_model(p)
    assert nesz.max() <= p.spec.nesz_db


def test_slick_detector_finds_strong_slick_and_ignores_open_sea():
    I, truth, la = slick_scene(contrast_db=8, seed=3)
    m = detect_dark(I)
    sea = ~truth & ~la
    assert (m & truth).sum() / truth.sum() > 0.8
    assert (m & sea).sum() / sea.sum() < 1e-3


def test_slick_detector_misses_weak_slick():
    I, truth, _ = slick_scene(contrast_db=2, seed=3)
    assert (detect_dark(I) & truth).sum() / truth.sum() < 0.05
