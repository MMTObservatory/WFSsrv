#!/usr/bin/env python

"""
MMT WFS Server
"""

import io
import os
import json
import pathlib
import redis

import numpy as np

import astropy.units as u

import logging
import logging.handlers

try:
    import tornado
except ImportError:
    raise RuntimeError("This server requires tornado.")
import tornado.web
import tornado.httpserver
import tornado.gen
import tornado.ioloop
import tornado.websocket
from tornado.process import Subprocess
from tornado.concurrent import run_on_executor
from tornado.log import enable_pretty_logging
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use('webagg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_webagg_core import (FigureCanvasWebAggCore, new_figure_manager_given_figure)

from mmtwfs.wfs import WFSFactory
from mmtwfs.zernike import ZernikeVector
from mmtwfs.telescope import MMT


glog = logging.getLogger('')
log = logging.getLogger('WFSsrv')
log.setLevel(logging.DEBUG)


def create_default_figures():
    zv = ZernikeVector(Z04=1)
    figures = {}
    ax = {}
    data = np.zeros((512, 512))
    tel = MMT(secondary='f5')  # secondary doesn't matter, just need methods for mirror forces/plots
    forces = tel.bending_forces(zv=zv)

    # stub for plot showing bkg-subtracted WFS image with aperture positions
    figures['slopes'], ax['slopes'] = plt.subplots()
    figures['slopes'].set_label("Aperture Positions and Spot Movement")
    ax['slopes'].imshow(data, cmap='Greys', origin='lower', interpolation='None')

    # stub for plot showing bkg-subtracted WFS image and residuals slopes of wavefront fit
    figures['residuals'], ax['residuals'] = plt.subplots()
    figures['residuals'].set_label("Zernike Fit Residuals")
    ax['residuals'].imshow(data, cmap='Greys', origin='lower', interpolation='None')

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
            self.render("home.html", current=self.application.wfs, wfslist=self.application.wfs_names)

    class SelectHandler(tornado.web.RequestHandler):
        def get(self):
            try:
                wfs = self.get_argument("wfs")
                if wfs in self.application.wfs_keys:
                    log.info(f"Setting {wfs}")
                    self.application.wfs = self.application.wfs_systems[wfs]
            except:
                log.warning(f"Must specify valid wfs: {wfs}.")
            finally:
                self.finish()

    class WFSPageHandler(tornado.web.RequestHandler):
        def get(self):
            try:
                wfs = self.get_argument("wfs")
                if wfs in self.application.wfs_keys:
                    log.info(f"Setting {wfs}")
                    self.application.wfs = self.application.wfs_systems[wfs]
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
                        "wfs.html",
                        wfsname=self.application.wfs.name,
                        ws_uris=ws_uris,
                        fig_ids=fig_ids,
                        figures=figkeys,
                        datadir=str(self.application.datadir) + "/",
                        modes=self.application.wfs.modes,
                        default_mode=self.application.wfs.default_mode,
                        m1_gain=self.application.wfs.m1_gain,
                        m2_gain=self.application.wfs.m2_gain,
                        log_uri=log_uri
                    )
            except Exception as e:
                log.warning(f"Must specify valid wfs: {wfs}. ({e})")
                self.finish()

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
        executor = ThreadPoolExecutor(max_workers=10)

        @tornado.gen.coroutine
        def async_plot(self, func, *args):
            k = yield func(*args)
            self.complete_refresh(k)

        @run_on_executor
        def make_psf(self, zernikes, telescope):
            log.debug("Making PSF image...")
            psf, self.application.figures['psf'] = telescope.psf(zv=zernikes.copy())
            return 'psf'

        @run_on_executor
        def make_wfmap(self, zernikes):
            log.debug("Making wavefront map...")
            self.application.figures['wavefront'] = zernikes.plot_map()
            return 'wavefront'

        @run_on_executor
        def make_barchart(self, zernikes, zrms, residual):
            log.debug("Making bar chart...")
            rms_asec = zrms.value / self.application.wfs.tiltfactor * u.arcsec
            self.application.figures['barchart'] = zernikes.bar_chart(
                last_mode=21,
                residual=residual,
                title=f"Total Wavefront RMS: {zrms.round(1)} ({rms_asec.round(2)})"
            )
            return 'barchart'

        @run_on_executor
        def make_fringebarchart(self, zernikes, focus, cc_x, cc_y):
            log.debug("Making fringe bar chart...")
            self.application.figures['fringebarchart'] = zernikes.fringe_bar_chart(
                title="Focus: {0:0.1f}  CC_X: {1:0.1f}  CC_Y: {2:0.1f}".format(
                    focus,
                    cc_x,
                    cc_y,
                ),
                max_c=1500*u.nm,
                last_mode=21
            )
            return 'fringebarchart'

        @run_on_executor
        def make_totalforces(self, telescope, forces, m1focus):
            log.debug("Making total forces plot...")
            self.application.figures['totalforces'] = telescope.plot_forces(forces, m1focus)
            self.application.figures['totalforces'].set_label("Total M1 Actuator Forces")
            return 'totalforces'

        @run_on_executor
        def make_pendingforces(self, telescope, forces, m1focus, limit):
            log.debug("Making pending forces plot...")
            self.application.figures['forces'] = telescope.plot_forces(
                forces,
                m1focus,
                limit=limit
            )
            self.application.figures['forces'].set_label("Requested M1 Actuator Forces")
            return 'forces'

        def complete_refresh(self, key):
            self.application.refresh_figure(key, self.application.figures[key])

        def get(self):
            self.application.close_figures()
            try:
                filename = self.get_argument("fitsfile")
                log.info(f"Analyzing {filename}...")
            except:
                log.warning("no wfs or file specified.")

            mode = self.get_argument("mode", default=None)
            connect = self.get_argument("connect", default=True)
            spher = self.get_argument("spher", default=False)

            if spher == "true":
                spher_mask = ['Z11', 'Z22']
                log.info(f"Ignoring all spherical terms {spher_mask}...")
            else:
                spher_mask = ['Z22']
                log.info(f"Only ignoring the high-order spherical terms {spher_mask}...")

            if os.path.isfile(filename) and not self.application.busy:
                self.application.busy = True
                if connect == "true":
                    self.application.wfs.connect()
                else:
                    self.application.wfs.disconnect()

                log.debug("Measuring slopes...")
                results = self.application.wfs.measure_slopes(filename, mode=mode, plot=True)
                if results['slopes'] is not None:
                    if 'seeing' in results:
                        log.info(f"Seeing (zenith): {results['seeing'].round(2)}")
                        log.info(f"Seeing (raw): {results['raw_seeing'].round(2)}")
                        if self.application.wfs.connected:
                            log.info("Publishing seeing values to redis.")
                            self.application.update_seeing(results)
                    tel = self.application.wfs.telescope

                    log.debug("Making slopes plot...")
                    self.application.figures['slopes'] = results['figures']['slopes']
                    self.application.refresh_figure('slopes', self.application.figures['slopes'])

                    log.debug("Fitting wavefront...")
                    zresults = self.application.wfs.fit_wavefront(results, plot=True)
                    log.info(f"Residual RMS: {zresults['residual_rms'].round(2)}")
                    self.application.figures['residuals'] = zresults['resid_plot']
                    self.application.refresh_figure('residuals', self.application.figures['residuals'])

                    zvec = zresults['zernike']
                    zvec_raw = zresults['rot_zernike']
                    zvec_ref = zresults['ref_zernike']

                    self.async_plot(self.make_psf, zvec.copy(), tel)
                    self.async_plot(self.make_wfmap, zvec.copy())

                    m1gain = self.application.wfs.m1_gain

                    # this is the total if we try to correct everything as fit
                    totforces, totm1focus, zv_totmasked = tel.calculate_primary_corrections(zvec, gain=m1gain)

                    self.async_plot(self.make_barchart, zvec.copy(), zresults['zernike_rms'], zresults['residual_rms'])
                    self.async_plot(self.make_totalforces, tel, totforces, totm1focus)

                    log.debug("Saving files and calculating corrections...")
                    zvec_file = self.application.datadir / (filename + ".zernike")
                    zvec_raw_file = self.application.datadir / (filename + ".raw.zernike")
                    zvec_ref_file = self.application.datadir / (filename + ".ref.zernike")
                    zvec.save(filename=zvec_file)
                    zvec_raw.save(filename=zvec_raw_file)
                    zvec_ref.save(filename=zvec_ref_file)
                    self.application.wavefront_fit = zvec

                    # check the RMS of the wavefront fit and only apply corrections if the fit is good enough.
                    # M2 can be more lenient to take care of large amounts of focus or coma.
                    if zresults['residual_rms'] < 450 * u.nm:
                        self.application.has_pending_m1 = True
                        self.application.has_pending_coma = True
                        self.application.has_pending_focus = True
                        log.info(f"{filename}: all proposed corrections valid.")
                    elif zresults['residual_rms'] <= 800 * u.nm:
                        self.application.has_pending_focus = True
                        log.warning(f"{filename}: only focus corrections valid.")
                    elif zresults['residual_rms'] > 800 * u.nm:
                        log.error(f"{filename}: wavefront fit too poor; no valid corrections")

                    self.application.has_pending_recenter = True

                    self.application.pending_focus = self.application.wfs.calculate_focus(zvec)

                    # only allow M1 corrections if we are reasonably close to good focus...
                    if self.application.pending_focus > 150 * u.um:
                        self.application.has_pending_m1 = False


                    self.application.pending_cc_x, self.application.pending_cc_y = self.application.wfs.calculate_cc(zvec)
                    self.async_plot(
                        self.make_fringebarchart,
                        zvec.copy(),
                        self.application.pending_focus,
                        self.application.pending_cc_x,
                        self.application.pending_cc_y
                    )

                    log.debug("Calculating pending forces...")
                    self.application.pending_az, self.application.pending_el = self.application.wfs.calculate_recenter(results)
                    self.application.pending_forces, self.application.pending_m1focus, zv_masked = \
                        self.application.wfs.calculate_primary(zvec, threshold=0.5*zresults['residual_rms'], mask=spher_mask)
                    self.application.pending_forcefile = self.application.datadir / (filename + ".forces")
                    zvec_masked_file = self.application.datadir / (filename + ".masked.zernike")
                    zv_masked.save(filename=zvec_masked_file)
                    limit = np.round(np.abs(self.application.pending_forces['force']).max())
                    self.async_plot(
                        self.make_pendingforces,
                        tel,
                        self.application.pending_forces,
                        self.application.pending_m1focus,
                        limit
                    )
                else:
                    log.error(f"Wavefront measurement failed: {filename}")
                    figures = create_default_figures()
                    figures['slopes'] = results['figures']['slopes']
                    self.application.refresh_figures(figures=figures)

            else:
                log.error(f"No such file: {filename}")

            self.application.wavefront_fit.denormalize()
            self.write(json.dumps(self.application.wavefront_fit.pretty_print()))
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
                fullfiles = sorted(p.glob("*_*.fits"), key=lambda x: x.stat().st_mtime)
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
            self.write(json.dumps(self.application.wavefront_fit.pretty_print()))
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

    class CompMirrorStatus(tornado.web.RequestHandler):
        def get(self):
            compmirror = self.application.wfs_systems['newf9'].compmirror
            status = compmirror.get_mirror()
            self.write(json.dumps(status))
            self.finish()

    class CompMirrorToggle(tornado.web.RequestHandler):
        def get(self):
            compmirror = self.application.wfs_systems['newf9'].compmirror
            status = compmirror.toggle_mirror()
            self.write(json.dumps(status))
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
            if "tornado.access" not in html and "poppy" not in html and "DEBUG" not in html:
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

    def restart_wfs(self, wfs):
        """
        If there's a configuration change, provide a way to reload the specified WFS
        """
        del self.wfs_systems[wfs]
        self.wfs_systems[wfs] = WFSFactory(wfs=wfs)

    def close_figures(self):
        if self.figures is not None:
            plt.close('all')

    def refresh_figure(self, k, figure):
        if k not in self.managers:
            fignum = id(figure)
            self.managers[k] = new_figure_manager_given_figure(fignum, figure)
            self.fig_id_map[fignum] = self.managers[k]
        else:
            canvas = FigureCanvasWebAggCore(figure)
            self.managers[k].canvas = canvas
            self.managers[k].canvas.manager = self.managers[k]
            self.managers[k]._get_toolbar(canvas)
            self.managers[k]._send_event("refresh")
            self.managers[k].canvas.draw()

    def refresh_figures(self, figures=None):
        if figures is None:
            self.figures = create_default_figures()
        else:
            self.figures = figures

        for k, f in self.figures.items():
            self.refresh_figure(k, f)

    def set_redis(self, key, value):
        """
        Set and publish redis 'key' to 'value'
        """
        resp = (None, None)
        try:
            r1 = self.redis_server.set(key, value)
            r2 = self.redis_server.publish(key, value)
            resp = (r1, r2)
        except Exception as e:
            log.warning(f"Failed to set redis {key} to {value}: {e}")
        return resp

    def update_seeing(self, results):
        try:
            wfs_seeing = results['seeing'].round(2).value
            wfs_raw_seeing = results['raw_seeing'].round(2).value
            r1 = self.set_redis('wfs_seeing', wfs_seeing)
            r2 = self.set_redis('wfs_raw_seeing', wfs_raw_seeing)
            if None not in r1 and None not in r2:
                log.info(f"Set redis values wfs_seeing={wfs_seeing} and wfs_raw_seeing={wfs_raw_seeing}")
            else:
                log.warning("Problem sending seeing values to redis...")
        except Exception as e:
            log.warning(f'Error connecting to MMTO API server... : {e}')

    def __init__(self):
        if 'WFSROOT' in os.environ:
            self.datadir = pathlib.Path(os.environ['WFSROOT'])
        elif 'HOME' in os.environ:
            self.datadir = pathlib.Path(os.environ['HOME']) / "wfsdat"
        else:
            self.datadir = pathlib.Path("wfsdat")

        if not self.datadir.exists():
            self.datadir.mkdir(parents=True)

        if not self.datadir.is_dir():
            self.datadir.unlink()
            self.datadir.mkdir(parents=True)

        if os.path.isdir(self.datadir):
            self.logfile = self.datadir / "wfs.log"
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            handler = logging.handlers.WatchedFileHandler(self.logfile)
            handler.setFormatter(formatter)
            glog.addHandler(handler)
            enable_pretty_logging()
        else:
            self.logfile = pathlib.Path("/dev/null")

        self.wfs = None
        self.wfs_systems = {}
        self.wfs_keys = ['newf9', 'f5', 'mmirs', 'binospec']
        self.wfs_names = {}
        for w in self.wfs_keys:
            self.wfs_systems[w] = WFSFactory(wfs=w)
            self.wfs_names[w] = self.wfs_systems[w].name

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

        self.redis_server = redis.StrictRedis(host='localhost', port=6379, db=0)

        handlers = [
            (r"/", self.HomeHandler),
            (r"/mpl\.js", tornado.web.RedirectHandler, dict(url="static/js/mpl.js")),
            (r"/select", self.SelectHandler),
            (r"/wfspage", self.WFSPageHandler),
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
            (r"/clear", self.ClearHandler),
            (r"/clearm1", self.ClearM1Handler),
            (r"/clearm2", self.ClearM2Handler),
            (r"/compmirror", self.CompMirrorStatus),
            (r"/compmirrortoggle", self.CompMirrorToggle),
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
    http_server.listen(8080)

    print("http://127.0.0.1:8080/")
    print("Press Ctrl+C to quit")

    tornado.ioloop.IOLoop.instance().start()
