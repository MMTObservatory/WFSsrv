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

**Key dependencies** (not on PyPI — installed from GitHub):
- `mmtwfs` — core wavefront analysis (slope measurement, Zernike fitting, mirror correction math)
- `camsrv` — camera server integration

**Environment variables**:
- `WFSROOT` — root directory for analysis output and logs (defaults to cwd)
- `REDISHOST` — Redis server for publishing results (default: `redis.mmto.arizona.edu`)

**Output file types** written to `WFSROOT`:
- `.zernike` — fitted Zernike coefficients
- `.forces` — M1 actuator force files
- `wfs.log` — application log

**Templates** (`wfssrv/templates/`): `home.html` (system selection), `wfs.html` (main analysis UI), `cwfs.html` (companion mirror UI for F/9).

**WebSocket handlers**: `LogStreamer` streams log output to the browser; `WebSocket` handles matplotlib figure interactivity.
