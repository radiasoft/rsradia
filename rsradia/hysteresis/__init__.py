# -*- coding: utf-8 -*-

u""":mod:`hysteresis` module of the rsradia package
"""

# Define number of interval divisions applied
# during Gauss-Kronrod quadrature integration
# (only used for anisotropic hysteresis)
KRON_SPLITS = 6

# Import implemented hysteresis models
from .jiles_atherton import JilesAtherton
from .preisach import Preisach