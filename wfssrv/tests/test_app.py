# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..wfssrv import WFSsrv as wfs


def test_wfs():
    app = wfs()
    assert app is not None
