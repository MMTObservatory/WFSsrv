# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..wfssrv import WFSsrv as wfs
from ..cwfssrv import WFSsrv as cwfs


def test_wfs():
    app = wfs()
    assert(app is not None)


def test_cwfs():
    app = cwfs()
    assert(app is not None)
