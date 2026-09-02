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
tox -e py313-test               # via tox
tox -e py313-cov                # with coverage (what CI runs)
```
The test suite is minimal — it just instantiates `WFSsrv()`. Because the constructor builds all four WFS systems via `mmtwfs`'s `WFSFactory`, even this smoke test requires `mmtwfs` and `camsrv` installed.

Python 3.13 is the minimum (`mmtwfs` requires >= 3.13); `tox.ini` defines `py313` and `py314` envs only. CI (`.github/workflows/wfssrv-tests.yml`) runs `py{313,314}-{cov,astropydev,numpydev}` plus `build_docs`, `linkcheck`, and `codestyle`.

The browser-side javascript in the templates is tested separately, with node's built-in test runner:
```bash
npm install   # once; jsdom is the only dependency
npm test      # runs wfssrv/tests/js/*.test.js
```
These tests read the inline `<script>` block out of `wfs.html` and `cwfs.html` and exercise it under jsdom, so they track what the server actually serves rather than a copy. The DOM fixtures in `wfssrv/tests/js/harness.js` are maintained by hand and guarded by a test asserting that every element the template script reaches for is present — if you add or rename a control, update the fixture. Node is only needed to run these tests; the served page has no build step. CI runs them from `.github/workflows/js-tests.yml` on node 22 and 24.

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

### Deployment

The summit installs from GitHub, into the `mmtwfs` conda environment used there:
```bash
pip install git+https://github.com/MMTObservatory/wfssrv#egg=wfssrv --upgrade
```
Then restart the server, and have the telescope operator reload their browser page to pick up template changes.

This installs from `master`, and the summit does not run from a working checkout — a fix is not deployed until it is pushed to `master`. Note that the templates ship as package data, so front-end changes reach the summit through this same install, not through a file copy.

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

**Templates** (`wfssrv/templates/`): `home.html` (system selection), `wfs.html` (main analysis UI), `cwfs.html` (comparison mirror UI for F/9).

**Front end**: Bootstrap 5.3.3 and its bundled JS are vendored under `wfssrv/static/` and served locally — there is no CDN dependency and **no jQuery**. All client code is plain ES: `fetch()` for the JSON/text endpoints, `addEventListener` for handlers, `data-bs-*` attributes for Bootstrap behaviour (tabs, modals). Keep it that way; do not reintroduce jQuery or `$()` idioms when editing the templates. (`plotly-latest.min.js` is still vendored and loaded by `wfs.html` but nothing calls it — it is a leftover.)

**WebSocket handlers**: `LogStreamer` streams log output to the browser (`/log`, consumed by a raw `new WebSocket(...)` in `wfs.html`/`cwfs.html`); `WebSocket` (`/<figure>/ws`) handles matplotlib figure interactivity.

**HTTP endpoints** (all registered in `WFSsrv.__init__`, all nested handler classes in `wfssrv/wfssrv.py`):
- Page/setup: `/`, `/select`, `/wfspage`, `/connect`, `/disconnect`, `/restart`, `/setdatadir`
- Analysis: `/analyze`, `/zfit`, `/files`, `/download_<type>.<ext>`
- Corrections: `/m1correct`, `/focuscorrect`, `/comacorrect`, `/recenter`, `/clear`, `/clearm1`, `/clearm2`, `/clearpending`
- Gains: `/m1gain`, `/m2gain`
- F/9 comparison mirror: `/compmirror`, `/compmirrortoggle`
- Matplotlib plumbing: `/mpl.js`, `/_static/*`, `/_images/*`

Most handlers accept `?format=json` and otherwise return plain text, which is what the templates rely on.

## Working with mmtwfs

`mmtwfs` is a git dependency, not a PyPI release, and is usually checked out alongside this repo (`../mmtwfs`). The APIs this server depends on are narrow: `WFSFactory`, `ZernikeVector` (`.bar_chart()`, `.fringe_bar_chart()`), `MMT` (`.bending_forces()`, `.calculate_primary_corrections()`, `.correct_primary()`, `.plot_forces()`), and `secondary.{focus,correct_coma,recenter}()`. When `mmtwfs` changes, those are the call sites to check — the smoke test in `wfssrv/tests/test_app.py` catches import/constructor breakage but nothing past it.
