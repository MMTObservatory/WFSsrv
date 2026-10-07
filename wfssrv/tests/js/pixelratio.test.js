"use strict";

// Moving the browser window between a regular and a retina display changes window.devicePixelRatio. mpl.js reads
// it once, when the figure is made, so without help the server keeps rendering at the old ratio and the plots come
// out half (or double) size until the page is reloaded.

const test = require("node:test");
const assert = require("node:assert");
const fs = require("fs");
const path = require("path");
const { watchPixelRatio } = require("../../static/js/pixelratio.js");

function fakeWindow(ratio) {
    const queries = [];
    const resizeListeners = [];
    return {
        devicePixelRatio: ratio,
        queries,
        addEventListener(type, fn) {
            if (type === "resize") resizeListeners.push(fn);
        },
        // the window is resized; pass a ratio to change it without the media query firing
        resize(newRatio) {
            if (newRatio !== undefined) this.devicePixelRatio = newRatio;
            resizeListeners.forEach((fn) => fn({}));
        },
        matchMedia(query) {
            const mq = {
                query,
                listeners: [],
                addEventListener(type, fn) {
                    this.listeners.push(fn);
                },
                removeEventListener(type, fn) {
                    this.listeners = this.listeners.filter((l) => l !== fn);
                },
            };
            queries.push(mq);
            return mq;
        },
        // move the window to a display with a different pixel ratio
        moveTo(newRatio) {
            const current = queries.filter((q) => q.listeners.length > 0);
            this.devicePixelRatio = newRatio;
            current.forEach((q) => q.listeners.slice().forEach((fn) => fn({ matches: false })));
        },
    };
}

function fakeFigure({ width = "640px", height = "480px", open = true } = {}) {
    const attrs = {};
    return {
        ratio: 1,
        ws: { readyState: open ? 1 : 0 },
        canvas: { style: { width, height }, setAttribute: (k, v) => (attrs[k] = v), attrs },
        messages: [],
        resizes: [],
        send_message(type, props) {
            this.messages.push({ type, ...props });
        },
        request_resize(w, h) {
            this.resizes.push([w, h]);
        },
    };
}

test("a move to a retina display re-renders the plots at the new pixel ratio", () => {
    const win = fakeWindow(1);
    const fig = fakeFigure();
    watchPixelRatio([fig], win);
    assert.equal(win.queries[0].query, "(resolution: 1dppx)");

    win.moveTo(2);
    assert.equal(fig.ratio, 2);
    assert.deepEqual(fig.messages, [{ type: "set_device_pixel_ratio", device_pixel_ratio: 2 }]);
    assert.equal(fig.canvas.attrs.width, 1280);
    assert.equal(fig.canvas.attrs.height, 960);
    assert.deepEqual(fig.resizes, [[640, 480]]);
});

test("moving back again is noticed too", () => {
    const win = fakeWindow(1);
    const fig = fakeFigure();
    watchPixelRatio([fig], win);
    win.moveTo(2);
    win.moveTo(1);
    assert.equal(fig.ratio, 1);
    assert.equal(fig.messages.at(-1).device_pixel_ratio, 1);
    assert.equal(win.queries.at(-1).query, "(resolution: 1dppx)");
});

test("a ratio change seen only as a window resize is handled too", () => {
    // not every browser fires the media query change when the window moves to another display
    const win = fakeWindow(1);
    const fig = fakeFigure();
    watchPixelRatio([fig], win);
    win.resize(2);
    assert.equal(fig.ratio, 2);
    assert.equal(fig.messages.length, 1);
});

test("a move that fires both the media query and a resize updates once, and plain resizes do nothing", () => {
    const win = fakeWindow(1);
    const fig = fakeFigure();
    watchPixelRatio([fig], win);
    win.moveTo(2);
    win.resize();
    win.resize();
    assert.equal(fig.messages.length, 1);
    assert.equal(fig.resizes.length, 1);
});

test("hidden and disconnected figures are not resized", () => {
    const win = fakeWindow(1);
    const hidden = fakeFigure({ width: "0px", height: "0px" });
    const closed = fakeFigure({ open: false });
    watchPixelRatio([hidden, closed], win);
    win.moveTo(2);
    // a hidden tab's figure still learns the new ratio, and gets its size when the tab is shown
    assert.equal(hidden.messages.length, 1);
    assert.deepEqual(hidden.resizes, []);
    // a closed websocket can't be sent anything
    assert.deepEqual(closed.messages, []);
    assert.deepEqual(closed.resizes, []);
});

test("both templates load the watcher and hand it their figures", () => {
    for (const t of ["wfs.html", "cwfs.html"]) {
        const src = fs.readFileSync(path.join(__dirname, "..", "..", "templates", t), "utf8");
        assert.match(src, /static_url\("js\/pixelratio\.js"\)/, `${t} should load pixelratio.js`);
        assert.match(src, /figures\.push\(fig\)/, `${t} should collect its figures`);
        assert.match(src, /watchPixelRatio\(figures\)/, `${t} should watch its figures`);
    }
});
