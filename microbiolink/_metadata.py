#!/usr/bin/env python

#
# This file is part of the `microbiolink` Python module
#
# Copyright 2026
# Korcsmaros Group, Imperial College London
#
# File author(s): Leila Potari-Gul (l.potari-gul@imperial.ac.uk)
#
# Distributed under the BSD-2-Clause license
# See the file `LICENSE` or read a copy at
# https://opensource.org/license/bsd-2-clause
#

"""Package metadata (version, authors, etc)."""

__all__ = ['__version__', '__author__', '__license__']

import importlib.metadata

_FALLBACK_VERSION = '0.0.1'

try:
    __version__ = importlib.metadata.version('microbiolink')
except importlib.metadata.PackageNotFoundError:
    # Package not installed (e.g. running from source checkout)
    __version__ = _FALLBACK_VERSION

__author__ = 'Leila Potari-Gul'
__license__ = 'BSD-2-Clause'
