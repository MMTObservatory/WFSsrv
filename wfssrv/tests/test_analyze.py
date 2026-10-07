# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib.resources
import logging
import os
import pathlib
import shutil
import tempfile
from unittest.mock import patch

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from tornado.testing import AsyncHTTPTestCase

from mmtwfs.zernike import ZernikeVector

from ..wfssrv import WFSsrv

FITSFILE = importlib.resources.files("mmtwfs") / "data" / "test_data" / "mmirs_wfs_0150.fits"


def _focus_only_results():
    slopes_fig = plt.figure()
    slopes_fig.set_label("WFS Image")
    period_fig = plt.figure()
    period_fig.set_label("Grid Periodicity")
    return {
        "slopes": None,
        "mode": "mmirs",
        "focus_only": True,
        "method": "periodicity",
        "grid": {"scale": 0.995, "scale_err": 0.003, "snr": np.array([1310.0, 1400.0])},
        "zernike": ZernikeVector(Z04=-390.0),
        "pending_focus": 2.5 * u.um,
        "figures": {"slopes": slopes_fig, "periodicity": period_fig},
    }


class TestFocusOnly(AsyncHTTPTestCase):
    def get_app(self):
        self.datadir = tempfile.TemporaryDirectory()
        self.root_handlers = list(logging.getLogger("").handlers)
        with patch.dict(os.environ, {"WFSROOT": self.datadir.name}):
            app = WFSsrv()
        self.fitsfile = pathlib.Path(self.datadir.name) / FITSFILE.name
        shutil.copy(FITSFILE, self.fitsfile)
        app.wfs = app.wfs_systems["mmirs"]
        return app

    def tearDown(self):
        super().tearDown()
        plt.close("all")
        # the app logs to <datadir>/wfs.log via the root logger. drop that handler before the directory goes away,
        # or later tests log into a missing file.
        root = logging.getLogger("")
        for h in set(root.handlers) - set(self.root_handlers):
            root.removeHandler(h)
            h.close()
        self.datadir.cleanup()

    def test_focus_only_offers_focus_alone(self):
        app = self._app
        # left over from an earlier, fully analyzed image
        app.has_pending_m1 = app.has_pending_coma = app.has_pending_recenter = True
        with patch.object(app.wfs, "measure_slopes", return_value=_focus_only_results()):
            resp = self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert resp.code == 200
        assert app.has_pending_focus and app.pending_focus == 2.5 * u.um
        assert not (app.has_pending_m1 or app.has_pending_coma or app.has_pending_recenter)
        assert app.wavefront_fit["Z04"] == -390.0 * u.nm
        assert pathlib.Path(str(self.fitsfile) + ".periodicity.zernike").exists()
        assert app.figures["slopes"].get_label() == "Grid Periodicity"
