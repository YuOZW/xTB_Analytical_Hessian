#!/usr/bin/env python3
"""Minimal post-install smoke test for a locally built wheel."""
import numpy as np
from tblite.interface import Calculator

numbers = np.array([8, 1, 1])
positions = np.array([
    [0.0000000000, 0.0000000000, 0.0],
    [1.4323367300, 1.1071526600, 0.0],
    [-1.4323367300, 1.1071526600, 0.0],
])
calc = Calculator('GFN1-xTB', numbers, positions)
calc.set('temperature', 0.0)
res = calc.hessian()
h = res.get('hessian')
freq = res.get('frequencies')
assert h.shape == (9, 9), h.shape
assert np.isfinite(h).all()
assert np.isfinite(freq).all()
print('tblite wheel smoke test: OK')
print('Hessian shape:', h.shape)
print('Frequencies / cm^-1:', freq)
