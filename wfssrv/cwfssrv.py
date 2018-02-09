#!/usr/bin/env python

"""
MMT CWFS Server

To get curvature wfs data via binospec single-object guider do the following via MSG:

1 capture_curvature_images <exp time (s)> <foc offset (um)>

1 sub curvature_files
"""

import sys
import io
import os
import json
import pathlib

import numpy as np

import astropy.units as u
from astropy.io import fits

import logging
import logging.handlers

try:
    import tornado
except ImportError:
    raise RuntimeError("This server requires tornado.")
import tornado.web
import tornado.httpserver
import tornado.ioloop
import tornado.websocket
from tornado.process import Subprocess
from tornado.log import enable_pretty_logging

import matplotlib
matplotlib.use('webagg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_webagg_core import (FigureCanvasWebAggCore, new_figure_manager_given_figure)

from mmtwfs.wfs import WFSFactory, center_pupil
from mmtwfs.zernike import ZernikeVector
from mmtwfs.telescope import MMT

sys.path.append("/mmt/cwfs/python")
sys.path.append("/Users/tim/src/cwfs/python")
from lsst.cwfs.instrument import Instrument
from lsst.cwfs.algorithm import Algorithm
from lsst.cwfs.image import Image, readFile
from lsst.cwfs.tools import outParam, outZer4Up

glog = logging.getLogger('')
log = logging.getLogger('CWFSsrv')
log.setLevel(logging.INFO)


def create_default_figures():
    zv = ZernikeVector(Z04=1)
    figures = {}
    ax = {}
    data = np.zeros((256, 256))
    tel = MMT(secondary='f5')  # secondary doesn't matter here, just need methods for mirror forces/plots
    forces = tel.bending_forces(zv=zv)

    # stub for plot showing bkg-subtracted WFS image with aperture positions
    figures['intra'], ax['intra'] = plt.subplots()
    figures['intra'].set_label("Intra-focal Image")
    ax['intra'].imshow(data, cmap='Greys', origin='lower', interpolation='None')

    # stub for plot showing bkg-subtracted WFS image and residuals slopes of wavefront fit
    figures['extra'], ax['extra'] = plt.subplots()
    figures['extra'].set_label("Extra-focal Image")
    ax['extra'].imshow(data, cmap='Greys', origin='lower', interpolation='None')

    # stub for zernike wavefront map
    figures['wavefront'] = zv.plot_map()

    # stub for zernike bar chart
    figures['barchart'] = zv.bar_chart()

    # stub for zernike fringe bar chart
    figures['fringebarchart'] = zv.fringe_bar_chart()

    # stubs for mirror forces
    figures['forces'] = tel.plot_forces(forces)
    figures['forces'].set_label("Requested M1 Actuator Forces")

    # stubs for mirror forces
    figures['totalforces'] = tel.plot_forces(forces)
    figures['totalforces'].set_label("Total M1 Actuator Forces")

    # stub for psf
    psf, figures['psf'] = tel.psf(zv=zv)

    return figures


class WFSsrv(tornado.web.Application):
    class HomeHandler(tornado.web.RequestHandler):
        """
        Serves the main HTML page.
        """
        def get(self):
            figkeys = []
            ws_uris = []
            fig_ids = []
            log_uri = "ws://{req.host}/log".format(req=self.request)
            for k, f in self.application.figures.items():
                manager = self.application.managers[k]
                fig_ids.append(manager.num)
                figkeys.append(k)
                ws_uri = "ws://{req.host}/{figdiv}/ws".format(req=self.request, figdiv=k)
                ws_uris.append(ws_uri)

            self.render(
                "cwfs.html",
                wfsname=self.application.wfs.name,
                ws_uris=ws_uris,
                fig_ids=fig_ids,
                figures=figkeys,
                datadir=str(self.application.datadir) + "/",
                m1_gain=self.application.wfs.m1_gain,
                m2_gain=self.application.wfs.m2_gain,
                log_uri=log_uri
            )

    class ConnectHandler(tornado.web.RequestHandler):
        def get(self):
            if self.application.wfs is not None:
                if not self.application.wfs.connected:
                    self.application.wfs.connect()
                    if self.application.wfs.connected:
                        log.info(f"{self.application.wfs.name} is connected.")
                    else:
                        log.warning(f"Couldn't connect to {self.application.wfs.name}. Offline?")
                else:
                    log.info(f"{self.application.wfs.name} already connected")
            self.finish()

    class DisconnectHandler(tornado.web.RequestHandler):
        def get(self):
            if self.application.wfs is not None:
                if self.application.wfs.connected:
                    self.application.wfs.disconnect()
                    log.info(f"{self.application.wfs.name} is disconnected.")
                else:
                    log.info(f"{self.application.wfs.name} already disconnected")
            self.finish()

    class AnalyzeHandler(tornado.web.RequestHandler):
        def get(self):
            self.application.close_figures()
            try:
                filename1 = self.get_argument("fitsfile1")
                filename2 = self.get_argument("fitsfile2")
                log.info(f"Analyzing {filename1} and {filename2}...")
            except:
                log.warning("No or not enough CWFS files specified.")

            connect = self.get_argument("connect", default=True)
            spher = self.get_argument("spher", default=False)
            model = self.get_argument("model", default="onAxis")
            thresh = self.get_argument("thresh", default=150.0)
            thresh = float(thresh) * u.nm
            focoff = self.get_argument("focoff", default=1000.0)
            focoff = float(focoff)

            tel = self.application.wfs.telescope

            if spher == "true":
                spher_mask = ['Z11', 'Z22']
                log.info(f"Ignoring the spherical terms {spher_mask}...")
            else:
                spher_mask = []

            if os.path.isfile(filename1) and os.path.isfile(filename2) and not self.application.busy:
                self.application.busy = True
                if connect == "true":
                    self.application.wfs.connect()
                else:
                    self.application.wfs.disconnect()

                # get rotator and focus values from the headers, if available
                rots = []
                focusvals = []
                images = [filename1, filename2]
                arrays = []
                for image in images:
                    h = fits.open(image)
                    hdr = h[-1].header
                    if 'ROT' in hdr:
                        rots.append(hdr['ROT'])
                        rot = hdr['ROT']
                    else:
                        rot = 0.0
                    if 'FOCUS' in hdr:
                        focusvals.append(hdr['FOCUS'])
                    data = h[-1].data
                    orig_shape = data.shape
                    pup_mask = tel.pupil_mask(rotation=rot, size=120)
                    x, y, fig = center_pupil(data, pup_mask, threshold=0.8, plot=False)
                    cenx = int(np.round(x))
                    ceny = int(np.round(y))
                    data = data[ceny-96:ceny+96, cenx-96:cenx+96]
                    if len(h) > 1 or orig_shape == (516, 532):  # check if we're raw binospec SOG images
                        log.info(f"Found raw Binospec SOG image. Flipping {image} up/down.")
                        data = np.flipud(data)
                    arrays.append(data)

                if len(rots) > 0:
                    rot = np.array(rots).mean() * u.deg
                    log.info(f"Using rotator angle of {rot.round(2)}.")
                else:
                    log.warning("WARNING: No rotator information in headers. Assuming rotator angle of 0.0 degrees.")
                    rot = 0.0 * u.deg

                if len(focusvals) == 2:
                    focusvals = np.array(focusvals)
                    I1 = Image(arrays[np.argmin(focusvals)], [0, 0], Image.INTRA)
                    I2 = Image(arrays[np.argmax(focusvals)], [0, 0], Image.EXTRA)
                    intra_file = pathlib.Path(images[np.argmin(focusvals)]).name
                    extra_file = pathlib.Path(images[np.argmax(focusvals)]).name
                    focoff = focusvals.max() - focusvals.mean()
                    log.info(f"Using an M2 focus offset of +/- {focoff} um.")
                    log.info(f"Intra-focal image: {images[np.argmin(focusvals)]}")
                    log.info(f"Extra-focal image: {images[np.argmax(focusvals)]}")
                else:
                    I1 = Image(readFile(images[0]), [0, 0], Image.INTRA)
                    I2 = Image(readFile(images[1]), [0, 0], Image.EXTRA)
                    intra_file = pathlib.Path(images[0]).name
                    extra_file = pathlib.Path(images[1]).name
                    log.warning(f"WARNING: No focus information in image headers. Assuming M2 focus offset of +/- {focoff} um.")
                    log.warning(f"WARNING: Assuming intra-focal image is {image[0]}")
                    log.warning(f"WARNING: Assuming extra-focal image is {image[1]}")

                # use second image filename as the root for output files
                output = filename2

                # The pupil rotation in the single-object guider on binospec was determined to be 0 deg.
                rotation = 0 * u.deg - rot
                log.info(f"Total pupil rotation: {rotation.round(2)}")

                # load instrument and algorithm parameters
                inst = Instrument('mmto', I1.sizeinPix)

                # this is a MMTO hack. 0.0 doesn't work, but this will yield an annular zernike solution that is very close
                # to circular. the MMTO wfs code currently doesn't support annular zernikes for calculating corrections.
                inst.obscuration = 0.01

                # convert M2 focus offset in microns to meters of focus shift at the instrument focal plane
                inst.offset = focoff * 1.0e-6 * 18.8

                # set up fitting algorithm
                algo = Algorithm("exp", inst, 0)

                log.info("Fitting wavefront using 'exp' algorithm and the 'onAxis' model.")

                # run it
                algo.runIt(inst, I1, I2, 'onAxis')

                log.info("Wavefront fit complete.")

                # output parameters
                outParam(output + ".param", algo, inst, I1, I2, 'onAxis')

                # set up ZernikeVector, denormalize it, and then derotate it
                zvec = ZernikeVector()
                zvec.from_array(algo.zer4UpNm, modestart=4, normalized=True)
                zvec.denormalize()
                zvec.rotate(angle=-rotation)

                log.debug("\n" + repr(zvec))

                # output Zernikes 4 and up
                outZer4Up(algo.zer4UpNm, 'nm', output + ".raw.lsst.zernikes")
                zvec.save(filename=output + ".rot.zernikes")

                m1gain = self.application.wfs.m1_gain

                # this is the total if we try to correct everything as fit
                totforces, totm1focus, zv_totmasked = tel.calculate_primary_corrections(zvec, gain=m1gain)
                figures = {}
                ax = {}

                # show intra-focal image
                figures['intra'], ax['intra'] = plt.subplots()
                figures['intra'].set_label("Intra-focal Image")
                im1 = ax['intra'].imshow(I1.image, cmap='Greys', origin='lower', interpolation='None')
                ax['intra'].set_title(intra_file)
                cbar1 = figures['intra'].colorbar(im1)

                # show extra-focal image
                figures['extra'], ax['extra'] = plt.subplots()
                figures['extra'].set_label("Extra-focal Image")
                im2 = ax['extra'].imshow(I2.image, cmap='Greys', origin='lower', interpolation='None')
                ax['extra'].set_title(extra_file)
                cbar2 = figures['extra'].colorbar(im2)

                # show wavefront map
                figures['wavefront'], ax['wavefront'] = plt.subplots()
                wfim = ax['wavefront'].imshow(algo.Wconverge * 1.0e9, origin="lower", cmap='RdBu')
                cbar_wf = figures['wavefront'].colorbar(wfim)
                cbar_wf.set_label(zvec.units.name, rotation=0)

                # show RMS bar chart
                zernike_rms = zvec.rms
                rms_asec = zernike_rms.value / self.application.wfs.tiltfactor * u.arcsec
                log.info(f"Total wavefront RMS: {zernike_rms.round(2)}")
                figures['barchart'] = zvec.bar_chart(
                    title=f"Total Wavefront RMS: {zernike_rms.round(1)} ({rms_asec.round(2)})"
                )

                # show total forces and PSF
                figures['totalforces'] = tel.plot_forces(totforces, totm1focus)
                figures['totalforces'].set_label("Total M1 Actuator Forces")
                psf, figures['psf'] = tel.psf(zv=zvec.copy())

                self.application.wavefront_fit = zvec
                self.application.pending_focus = self.application.wfs.calculate_focus(zvec)
                self.application.pending_cc_x, self.application.pending_cc_y = self.application.wfs.calculate_cc(zvec)

                self.application.pending_forces, self.application.pending_m1focus, zv_masked = \
                    self.application.wfs.calculate_primary(zvec, threshold=thresh, mask=spher_mask)

                self.application.pending_forcefile = output + ".forces"
                zvec_masked_file = output + ".masked.zernike"
                zv_masked.save(filename=zvec_masked_file)

                limit = np.round(np.abs(self.application.pending_forces['force']).max())
                figures['forces'] = tel.plot_forces(
                    self.application.pending_forces,
                    self.application.pending_m1focus,
                    limit=limit
                )
                figures['forces'].set_label("Requested M1 Actuator Forces")

                figures['fringebarchart'] = zvec.fringe_bar_chart(
                    title="Focus: {0:0.1f}  CC_X: {1:0.1f}  CC_Y: {2:0.1f}".format(
                        self.application.pending_focus,
                        self.application.pending_cc_x,
                        self.application.pending_cc_y,
                    )
                )
                self.application.has_pending_m1 = True
                self.application.has_pending_focus = True
                self.application.has_pending_coma = True

                self.application.refresh_figures(figures=figures)
            else:
                log.error(f"No such files: {filename1} {filename2}")

            self.application.wavefront_fit.denormalize()
            self.write(json.dumps(repr(self.application.wavefront_fit)))
            self.application.busy = False
            self.finish()

    class M1CorrectHandler(tornado.web.RequestHandler):
        def get(self):
            log.info("M1 corrections...")
            if self.application.has_pending_m1 and self.application.wfs.connected:
                self.application.wfs.telescope.correct_primary(
                    self.application.pending_forces,
                    self.application.pending_m1focus,
                    filename=self.application.pending_forcefile
                )
                max_f = self.application.pending_forces['force'].max()
                min_f = self.application.pending_forces['force'].min()
                log.info(f"Maximum force {round(max_f, 2)} N")
                log.info(f"Minimum force {round(min_f, 2)} N")
                log.info("Adjusting M1 focus by {0:0.1f}".format(self.application.pending_m1focus))
                self.application.has_pending_m1 = False
                self.write("Sending forces to cell and {0:0.1f} focus to secondary...".format(self.application.pending_m1focus))
            else:
                log.info("No M1 corrections sent")
                self.write("No M1 corrections sent")
            self.finish()

    class FocusCorrectHandler(tornado.web.RequestHandler):
        def get(self):
            log.info("M2 focus corrections...")
            if self.application.has_pending_focus and self.application.wfs.connected:
                self.application.wfs.secondary.focus(self.application.pending_focus)
                self.application.has_pending_focus = False
                log_str = "Sending {0:0.1f} focus to secondary...".format(
                    self.application.pending_focus
                )
                log.info(log_str)
                self.write(log_str)
            else:
                log.info("No focus corrections sent")
                self.write("No focus corrections sent")
            self.finish()

    class ComaCorrectHandler(tornado.web.RequestHandler):
        def get(self):
            log.info("M2 coma corrections...")
            if self.application.has_pending_coma and self.application.wfs.connected:
                self.application.wfs.secondary.correct_coma(self.application.pending_cc_x, self.application.pending_cc_y)
                self.application.has_pending_coma = False
                log_str = "Sending {0:0.1f}/{1:0.1f} CC_X/CC_Y to secondary...".format(
                    self.application.pending_cc_x,
                    self.application.pending_cc_y
                )
                log.info(log_str)
                self.write(log_str)
            else:
                log.info("No coma corrections sent")
                self.write("No coma corrections sent")
            self.finish()

    class RecenterHandler(tornado.web.RequestHandler):
        def get(self):
            log.info("Recentering...")
            if self.application.has_pending_recenter and self.application.wfs.connected:
                self.application.wfs.secondary.recenter(self.application.pending_az, self.application.pending_el)
                self.application.has_pending_recenter = False
                log_str = "Sending {0:0.1f}/{1:0.1f} of az/el to recenter...".format(
                    self.application.pending_az,
                    self.application.pending_el
                )
                log.info(log_str)
                self.write(log_str)
            else:
                log.info("No M2 recenter corrections sent")
                self.write("No M2 recenter corrections sent")
            self.finish()

    class RestartHandler(tornado.web.RequestHandler):
        def get(self):
            try:
                wfs = self.get_argument('wfs')
                self.application.restart_wfs(wfs)
                log.info(f"restarting {wfs}")
            except:
                log.info("no wfs specified")
            finally:
                self.finish()

    class DataDirHandler(tornado.web.RequestHandler):
        def get(self):
            try:
                datadir = self.get_argument("datadir")
                if os.path.isdir(datadir):
                    log.info(f"setting datadir to {datadir}")
                    self.application.datadir = pathlib.Path(datadir)
            except:
                log.info("no datadir specified")
            finally:
                self.finish()

    class M1GainHandler(tornado.web.RequestHandler):
        def get(self):
            if self.application.wfs is not None:
                self.write(f"{self.application.wfs.m1_gain}")
            else:
                self.write("no WFS selected.")

        def post(self):
            self.set_header("Content-Type", "text/plain")
            try:
                gain = float(self.get_body_argument('gain'))
                if self.application.wfs is not None:
                    if gain >= 0.0 and gain <= 1.0:
                        self.application.wfs.m1_gain = gain
                    else:
                        log.warning(f"Invalid M1 gain: {gain}")
                    log.info(f"Setting M1 gain to {gain}")
            except Exception as e:
                log.warning(f"No M1 gain specified: {e}")
                log.info(f"Body: {self.request.body}")

    class M2GainHandler(tornado.web.RequestHandler):
        def get(self):
            if self.application.wfs is not None:
                self.write(f"{self.application.wfs.m2_gain}")
            else:
                self.write("no WFS selected.")
            self.finish()

        def post(self):
            self.set_header("Content-Type", "text/plain")
            try:
                gain = float(self.get_body_argument('gain'))
                if self.application.wfs is not None:
                    if gain >= 0.0 and gain <= 1.0:
                        self.application.wfs.m2_gain = gain
                    else:
                        log.warning(f"Invalid M2 gain: {gain}")
                    log.info(f"Setting M2 gain to {gain}")
            except Exception as e:
                log.warning(f"No M2 gain specified: {e}")
                log.info(f"Body: {self.request.body}")
            finally:
                self.finish()

    class PendingHandler(tornado.web.RequestHandler):
        def get(self):
            self.write(f"M1: {self.application.has_pending_m1}")
            self.write(f"M2: {self.application.has_pending_m2}")
            self.write(f"recenter: {self.application.has_pending_recenter}")
            self.finish()

        def post(self):
            self.application.has_pending_m1 = False
            self.application.has_pending_m2 = False
            self.application.has_pending_recenter = False
            self.finish()

    class FilesHandler(tornado.web.RequestHandler):
        def get(self):
            p = self.application.datadir
            try:
                fullfiles = sorted(p.glob("sog*_*.fits"), key=lambda x: x.stat().st_mtime)
            except PermissionError as e:
                # started getting weird permission errors on hacksaw that looks like NFS race bug.
                # running 'ls' in the directory clears the error...
                log.warning(f"Permission error while listing files in {p}...")
                os.system(f"ls {p} > /dev/null")
                fullfiles = []
            files = []
            for f in fullfiles:
                files.append(f.name)
            files.reverse()
            self.write(json.dumps(files))
            self.finish()

    class ZernikeFitHandler(tornado.web.RequestHandler):
        def get(self):
            self.application.wavefront_fit.denormalize()
            self.write(json.dumps(repr(self.application.wavefront_fit)))
            self.finish()

    class ClearM1Handler(tornado.web.RequestHandler):
        def get(self):
            self.application.wfs.clear_m1_corrections()
            log_str = "Cleared M1 forces and M2 m1spherical offsets..."
            self.write(log_str)
            self.finish()

    class ClearM2Handler(tornado.web.RequestHandler):
        def get(self):
            self.application.wfs.clear_m2_corrections()
            log_str = "Cleared M2 wfs offsets..."
            self.write(log_str)
            self.finish()

    class ClearHandler(tornado.web.RequestHandler):
        def get(self):
            self.application.close_figures()
            self.application.wfs.clear_corrections()
            figures = create_default_figures()
            self.application.refresh_figures(figures=figures)
            log_str = "Cleared M1 forces and M2 wfs/m1spherical offsets...."
            self.write(log_str)
            self.finish()

    class Download(tornado.web.RequestHandler):
        """
        Handles downloading of the figure in various file formats.
        """
        def get(self, fig, fmt):
            managers = self.application.managers

            mimetypes = {
                'ps': 'application/postscript',
                'eps': 'application/postscript',
                'pdf': 'application/pdf',
                'svg': 'image/svg+xml',
                'png': 'image/png',
                'jpeg': 'image/jpeg',
                'tif': 'image/tiff',
                'emf': 'application/emf'
            }

            self.set_header('Content-Type', mimetypes.get(fmt, 'binary'))

            buff = io.BytesIO()
            managers[fig].canvas.print_figure(buff, format=fmt)
            self.write(buff.getvalue())
            self.finish()

    class LogStreamer(tornado.websocket.WebSocketHandler):
        """
        A websocket for streaming log messages from log file to the browser.
        """
        def open(self):
            filename = self.application.logfile
            self.proc = Subprocess(["tail", "-f", "-n", "0", filename],
                                   stdout=Subprocess.STREAM,
                                   bufsize=1)
            self.proc.set_exit_callback(self._close)
            self.proc.stdout.read_until(b"\n", self.write_line)

        def _close(self, *args, **kwargs):
            self.close()

        def on_close(self, *args, **kwargs):
            log.info("Trying to kill log streaming process...")
            self.proc.proc.terminate()
            self.proc.proc.wait()

        def write_line(self, data):
            html = data.decode()
            if "WARNING" in html:
                color = "text-warning"
            elif "ERROR" in html:
                color = "text-danger"
            else:
                color = "text-success"
            if "tornado.access" not in html and "poppy" not in html:
                html = "<samp><span class=%s>%s</span></samp>" % (color, html)
                html += "<script>$(\"#log\").scrollTop($(\"#log\")[0].scrollHeight);</script>"
                self.write_message(html.encode())
            self.proc.stdout.read_until(b"\n", self.write_line)

    class WebSocket(tornado.websocket.WebSocketHandler):
        """
        A websocket for interactive communication between the plot in
        the browser and the server.

        In addition to the methods required by tornado, it is required to
        have two callback methods:

            - ``send_json(json_content)`` is called by matplotlib when
              it needs to send json to the browser.  `json_content` is
              a JSON tree (Python dictionary), and it is the responsibility
              of this implementation to encode it as a string to send over
              the socket.

            - ``send_binary(blob)`` is called to send binary image data
              to the browser.
        """
        supports_binary = True

        def open(self, figname):
            # Register the websocket with the FigureManager.
            self.figname = figname
            manager = self.application.managers[figname]
            manager.add_web_socket(self)
            if hasattr(self, 'set_nodelay'):
                self.set_nodelay(True)

        def on_close(self):
            # When the socket is closed, deregister the websocket with
            # the FigureManager.
            manager = self.application.managers[self.figname]
            manager.remove_web_socket(self)

        def on_message(self, message):
            # The 'supports_binary' message is relevant to the
            # websocket itself.  The other messages get passed along
            # to matplotlib as-is.

            # Every message has a "type" and a "figure_id".
            message = json.loads(message)
            if message['type'] == 'supports_binary':
                self.supports_binary = message['value']
            else:
                manager = self.application.fig_id_map[message['figure_id']]
                manager.handle_json(message)

        def send_json(self, content):
            self.write_message(json.dumps(content))

        def send_binary(self, blob):
            if self.supports_binary:
                self.write_message(blob, binary=True)
            else:
                data_uri = "data:image/png;base64,{0}".format(
                    blob.encode('base64').replace('\n', ''))
                self.write_message(data_uri)

    def close_figures(self):
        if self.figures is not None:
            plt.close('all')

    def refresh_figures(self, figures=None):
        if figures is None:
            self.figures = create_default_figures()
        else:
            self.figures = figures

        for k, f in self.figures.items():
            if k not in self.managers:
                fignum = id(f)
                self.managers[k] = new_figure_manager_given_figure(fignum, f)
                self.fig_id_map[fignum] = self.managers[k]
            else:
                canvas = FigureCanvasWebAggCore(f)
                self.managers[k].canvas = canvas
                self.managers[k].canvas.manager = self.managers[k]
                self.managers[k]._get_toolbar(canvas)
                self.managers[k]._send_event("refresh")
                self.managers[k].canvas.draw()

    def __init__(self):
        if 'CWFSROOT' in os.environ:
            self.datadir = pathlib.Path(os.environ['CWFSROOT'])
        else:
            self.datadir = pathlib.Path("/mmt/shwfs/datadir")

        if os.path.isdir(self.datadir):
            self.logfile = self.datadir / "cwfs.log"
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            handler = logging.handlers.WatchedFileHandler(self.logfile)
            handler.setFormatter(formatter)
            glog.addHandler(handler)
            enable_pretty_logging()
        else:
            self.logfile = pathlib.Path("/dev/null")

        self.wfs = WFSFactory(wfs='binospec')  # CWFS currently is only used with binospec
        self.wfs.name = "CWFS"

        self.busy = False
        self.has_pending_m1 = False
        self.has_pending_coma = False
        self.has_pending_focus = False
        self.has_pending_recenter = False

        self.figures = None
        self.managers = {}
        self.fig_id_map = {}
        self.refresh_figures()
        self.wavefront_fit = ZernikeVector(Z04=1)

        handlers = [
            (r"/", self.HomeHandler),
            (r"/mpl\.js", tornado.web.RedirectHandler, dict(url="static/js/mpl.js")),
            (r"/connect", self.ConnectHandler),
            (r"/disconnect", self.DisconnectHandler),
            (r"/analyze", self.AnalyzeHandler),
            (r"/m1correct", self.M1CorrectHandler),
            (r"/focuscorrect", self.FocusCorrectHandler),
            (r"/comacorrect", self.ComaCorrectHandler),
            (r"/recenter", self.RecenterHandler),
            (r"/restart", self.RestartHandler),
            (r"/setdatadir", self.DataDirHandler),
            (r"/m1gain", self.M1GainHandler),
            (r"/m2gain", self.M2GainHandler),
            (r"/clearpending", self.PendingHandler),
            (r"/files", self.FilesHandler),
            (r"/zfit", self.ZernikeFitHandler),
            (r"/clearm1", self.ClearM1Handler),
            (r"/clearm2", self.ClearM2Handler),
            (r"/clear", self.ClearHandler),
            (r'/download_([a-z]+).([a-z0-9.]+)', self.Download),
            (r'/log', self.LogStreamer),
            (r'/([a-z0-9.]+)/ws', self.WebSocket)
        ]

        settings = dict(
            template_path=os.path.join(os.path.dirname(__file__), "templates"),
            static_path=os.path.join(os.path.dirname(__file__), "static"),
            debug=True
        )
        super(WFSsrv, self).__init__(handlers, **settings)


if __name__ == "__main__":
    application = WFSsrv()

    http_server = tornado.httpserver.HTTPServer(application)
    http_server.listen(8081)

    print("http://127.0.0.1:8081/")
    print("Press Ctrl+C to quit")

    tornado.ioloop.IOLoop.instance().start()
