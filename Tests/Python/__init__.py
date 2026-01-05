# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import sys
from contextlib import contextmanager


@contextmanager
def push_argv(argv):
    old_argv = sys.argv
    sys.argv = argv
    yield
    sys.argv = old_argv
