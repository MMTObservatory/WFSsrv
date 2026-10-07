# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging
import os
import tempfile
from unittest.mock import patch

import matplotlib.pyplot as plt
from tornado.testing import AsyncHTTPTestCase

from ..wfssrv import WFSsrv

FLAGS = ["has_pending_m1", "has_pending_focus", "has_pending_coma", "has_pending_recenter"]


class TestClearPending(AsyncHTTPTestCase):
    def get_app(self):
        plt.close("all")  # each app opens its own figures; don't inherit earlier tests'
        self.datadir = tempfile.TemporaryDirectory()
        self.root_handlers = list(logging.getLogger("").handlers)
        with patch.dict(os.environ, {"WFSROOT": self.datadir.name}):
            return WFSsrv()

    def tearDown(self):
        super().tearDown()
        plt.close("all")
        # drop the app's <datadir>/wfs.log handler before the directory goes away
        root = logging.getLogger("")
        for h in set(root.handlers) - set(self.root_handlers):
            root.removeHandler(h)
            h.close()
        self.datadir.cleanup()

    def test_get_reports_pending(self):
        self._app.has_pending_coma = True
        resp = self.fetch("/clearpending")
        assert resp.code == 200
        assert resp.body.decode().split("\n") == ["M1: False", "focus: False", "coma: True", "recenter: False"]

    def test_post_clears_pending(self):
        for f in FLAGS:
            setattr(self._app, f, True)
        resp = self.fetch("/clearpending", method="POST", body="")
        assert resp.code == 200
        assert not any(getattr(self._app, f) for f in FLAGS)
