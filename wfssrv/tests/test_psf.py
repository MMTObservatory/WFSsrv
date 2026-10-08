# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging
import os
import pathlib
import shutil
import tempfile
from unittest.mock import patch

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.text import Text
from tornado.testing import AsyncHTTPTestCase

from mmtwfs.psf import PSF_BANDS

from ..wfssrv import WFSsrv
from .test_analyze import FITSFILE, _focus_only_results


def _titles(fig):
    return " ".join(t.get_text() for t in fig.findobj(Text))


class TestPSF(AsyncHTTPTestCase):
    def get_app(self):
        plt.close("all")
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
        root = logging.getLogger("")
        for h in set(root.handlers) - set(self.root_handlers):
            root.removeHandler(h)
            h.close()
        self.datadir.cleanup()

    def test_page_has_psf_panel_and_band_menu(self):
        app = self._app
        assert "psf" in app.figures
        assert app.psf_band == "500nm"
        page = self.fetch("/wfspage?wfs=mmirs").body.decode()
        assert 'id="psf"' in page
        assert 'id="psfband"' in page
        for band in PSF_BANDS:
            assert f'value="{band}"' in page

    def test_full_analysis_keeps_wavefront_and_raw_seeing(self):
        app = self._app
        resp = self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert resp.code == 200
        assert app.psf_wavefront is not None
        assert app.psf_wavefront["Z04"] == app.wavefront_fit["Z04"]
        # the measured wavefront is as observed, so the delivered PSF uses the seeing at the observed airmass
        assert app.psf_seeing is not None
        assert u.isclose(app.psf_seeing, 0.926 * u.arcsec, atol=0.01 * u.arcsec)

        assert app.update_psf() == "psf"
        titles = _titles(app.figures["psf"])
        assert "PSF at 500nm" in titles
        assert "Optics" in titles and "Delivered" in titles

    def test_band_menu_redraws_psf(self):
        app = self._app
        self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        resp = self.fetch("/psfband?band=K")
        assert resp.code == 200
        assert app.psf_band == "K"
        assert "PSF at K" in _titles(app.figures["psf"])

    def test_bogus_band_rejected(self):
        app = self._app
        resp = self.fetch("/psfband?band=Q")
        assert resp.code == 400
        assert app.psf_band == "500nm"

    def test_band_choice_kept_without_wavefront(self):
        app = self._app
        stub = app.figures["psf"]
        resp = self.fetch("/psfband?band=J")
        assert resp.code == 200
        assert app.psf_band == "J"
        assert app.figures["psf"] is stub

    def test_focus_only_shows_optics_alone(self):
        app = self._app
        with patch.object(app.wfs, "measure_slopes", return_value=_focus_only_results()):
            self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert list(app.psf_wavefront.coeffs) == ["Z04"]
        assert app.psf_seeing is None
        app.update_psf()
        titles = _titles(app.figures["psf"])
        assert "Optics" in titles and "Delivered" not in titles

    def test_failed_analysis_forgets_old_wavefront(self):
        app = self._app
        self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert app.psf_wavefront is not None
        with patch.object(app.wfs, "measure_slopes", side_effect=ValueError("not a WFS image")):
            self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert app.psf_wavefront is None
        # changing the band must not bring back the PSF of the earlier image
        stub = app.figures["psf"]
        self.fetch("/psfband?band=H")
        assert app.figures["psf"] is stub

    def test_clear_forgets_wavefront(self):
        app = self._app
        self.fetch(f"/analyze?connect=false&fitsfile={self.fitsfile}")
        assert app.psf_wavefront is not None and app.psf_seeing is not None
        # stand in for the hexapod and cell, which clearing would otherwise command if connected
        with patch.object(app.wfs, "clear_corrections") as clear:
            resp = self.fetch("/clear")
        assert resp.code == 200
        clear.assert_called_once()
        assert app.psf_wavefront is None
        assert app.psf_seeing is None
        # the PSF panel is back to its placeholder and the band menu doesn't bring the old PSF back
        stub = app.figures["psf"]
        assert "Optics" not in _titles(stub)
        self.fetch("/psfband?band=H")
        assert app.figures["psf"] is stub
