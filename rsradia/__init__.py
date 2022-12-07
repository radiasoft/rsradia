# -*- coding: utf-8 -*-

u""":mod:`rsradia` package

:copyright: Copyright (c) 2020 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

from __future__ import absolute_import, division, print_function
import pkg_resources
from pykern import pkio
from pykern import pkresource

# Set the version number if one exists
try:
    # We only have a version once the package is installed.
    __version__ = pkg_resources.get_distribution('rsradia').version
except pkg_resources.DistributionNotFound:
    pass

# Define pi & the vacuum permeability constant
PI = 3.141592653589793
MU0 = 4*PI*1e-7

# TODO: Not sure how to handle dump paths. Currently just write and run in the notebook's directory.
_PARALLEL_RADIA_TEMPLATE = 'run_parallel_radia.py.jinja'
_PARALLEL_RADIA_SCRIPT = _PARALLEL_RADIA_TEMPLATE.rsplit('.', 1)[0]
_TEMPLATE_PATH = pkio.py_path(pkresource.filename(''))
_DEFAULT_SOLVER = 'RlxPre'
_NOTEBOOK_DUMP_PATH = './temp_mag.bin'
_SCRIPT_DUMP_PATH = './solver_result.bin'  # TODO: This has to go to both mpi_solve and the template
_MPI = 'mpiexec'