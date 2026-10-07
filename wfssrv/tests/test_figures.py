# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging
import os
import tempfile
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np

from ..wfssrv import WFSsrv


def test_swapped_figure_matches_browser_canvas():
    plt.close("all")
    datadir = tempfile.TemporaryDirectory()
    root_handlers = list(logging.getLogger("").handlers)
    try:
        with patch.dict(os.environ, {"WFSROOT": datadir.name}):
            app = WFSsrv()
        # the app already shows a figure in every panel; a retina browser connects and sizes it
        canvas = app.managers["slopes"].canvas
        canvas.figure.set_size_inches(5.0, 4.0, forward=False)
        canvas._set_device_pixel_ratio(2.0)
        shown = canvas.get_width_height(physical=True)

        # a new analysis hands over a figure of a different size. it has to come out the same size in pixels as
        # the one the browser is showing, or the browser draws a cropped, magnified image until it is resized.
        new = plt.figure(figsize=(8.0, 8.0))
        app.refresh_figure("slopes", new)
        canvas = app.managers["slopes"].canvas
        assert canvas.figure is new
        assert canvas.device_pixel_ratio == 2.0
        assert np.allclose(new.get_size_inches(), (5.0, 4.0))
        assert canvas.get_width_height(physical=True) == shown
    finally:
        root = logging.getLogger("")
        for h in set(root.handlers) - set(root_handlers):
            root.removeHandler(h)
            h.close()
        plt.close("all")
        datadir.cleanup()


def test_full_frame_after_browser_resize():
    # showing a hidden tab resizes (and so clears) the browser's canvas. if the figure is already that size, the next
    # frame used to be a diff against an image the browser no longer has, leaving the panel blank.
    plt.close("all")
    datadir = tempfile.TemporaryDirectory()
    root_handlers = list(logging.getLogger("").handlers)
    try:
        with patch.dict(os.environ, {"WFSROOT": datadir.name}):
            app = WFSsrv()
        canvas = app.managers["barchart"].canvas
        w, h = canvas.get_width_height()
        canvas.draw()
        canvas.get_diff_image()  # the frame the server sent while the tab was hidden
        canvas.handle_resize({"width": w, "height": h})
        canvas.draw()
        canvas.get_diff_image()
        assert canvas._current_image_mode == "full"
    finally:
        root = logging.getLogger("")
        for h in set(root.handlers) - set(root_handlers):
            root.removeHandler(h)
            h.close()
        plt.close("all")
        datadir.cleanup()
