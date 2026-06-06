# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WFSsrv is a Tornado-based web server for analyzing wavefront sensor (WFS) data from the MMT Observatory. It processes FITS images, fits Zernike polynomials, and provides controls for correcting the telescope's primary (M1) and secondary (M2) mirrors.

## Commands

### Running the server
```bash
wfssrv                  # installed entry point
pixi run wfssrv         # via pixi
```

### Testing
```bash
pytest                          # run all tests
pytest wfssrv/tests/test_app.py # run a single test file
pixi run test                   # via pixi
tox -e py312-test               # via tox
```
The test suite is minimal — it just instantiates `WFSsrv()`. Because the constructor builds all four WFS systems via `mmtwfs`'s `WFSFactory`, even this smoke test requires `mmtwfs` and `camsrv` installed.

### Code style
```bash
flake8 wfssrv --max-line-length=135   # lint (max line length is 135)
black wfssrv                           # format
```

### Documentation
```bash
tox -e build_docs   # build Sphinx docs
tox -e linkcheck    # check doc links
```

## Architecture

The entire application is a single `WFSsrv(tornado.web.Application)` class in `wfssrv/wfssrv.py`. All HTTP handlers are defined as nested classes within that file. The app listens on port 8080.

**Supported WFS systems**: `newf9`, `f5`, `mmirs`, `binospec` — selected at runtime via the web UI.

**Analysis pipeline** (triggered by `AnalyzeHandler`):
1. Reads a FITS image from disk
2. Delegates to `mmtwfs` library for slope measurement and Zernike fitting
3. Calculates M1 force corrections and M2 (focus/coma) adjustments
4. Streams result plots to the browser via matplotlib's WebAgg backend over WebSocket
5. Publishes seeing and correction data to Redis

**Correction state machine**: The app instance carries mutable state — `busy` (a single-flight lock; `AnalyzeHandler` ignores requests while set) and `has_pending_m1/coma/focus/recenter` flags. `AnalyzeHandler` gates which corrections become "pending" on the wavefront fit quality: residual RMS < 4000 nm enables all corrections, ≤ 7000 nm enables focus only, > 7000 nm enables none (M1 is additionally suppressed if pending focus > 150 µm). The `M1Correct`/`FocusCorrect`/`ComaCorrect`/`Recenter` handlers each apply their pending correction only if its flag is set *and* the WFS is connected, then clear the flag. **Hexapod recentering is currently disabled** — `RecenterHandler` only reports the mount offsets to apply manually.

**Figure streaming**: Six named figures (`slopes`, `residuals`, `barchart`, `fringebarchart`, `forces`, `totalforces`) are created up front by `create_default_figures()` and live in `self.figures`. Each has a persistent WebAgg `FigureManager` in `self.managers` (keyed by name) plus a `self.fig_id_map` (keyed by matplotlib figure id, used to route incoming WebSocket messages). `refresh_figure()` swaps a manager's canvas in place rather than recreating the manager, preserving the browser's WebSocket connection. The slower plots (bar charts, force plots) are rendered on a `ThreadPoolExecutor` and pushed to the browser asynchronously via `async_plot()`.

Each WFS system is hard-coded to fit `nzern = 10` Zernike modes.

**Key dependencies** (not on PyPI — installed from GitHub):
- `mmtwfs` — core wavefront analysis (slope measurement, Zernike fitting, mirror correction math)
- `camsrv` — camera server integration

**Environment variables**:
- `WFSROOT` — root directory for analysis output and logs (defaults to cwd)
- `REDISHOST` — Redis server for publishing results (default: `redis.mmto.arizona.edu`)

**Output file types** written to `WFSROOT` (the analysis `datadir`, settable at runtime via `/setdatadir`):
- `.zernike` — fitted Zernike coefficients (also `.raw.zernike`, `.ref.zernike`, `.masked.zernike` variants)
- `.forces` — M1 actuator force files
- `wfs.log` — application log (tailed live to the browser by `LogStreamer`)

**Templates** (`wfssrv/templates/`): `home.html` (system selection), `wfs.html` (main analysis UI), `cwfs.html` (companion mirror UI for F/9).

**WebSocket handlers**: `LogStreamer` streams log output to the browser; `WebSocket` handles matplotlib figure interactivity.
