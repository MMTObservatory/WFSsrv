"use strict";

// Test harness for the interactive javascript embedded in the wfssrv templates.
//
// The javascript under test is pulled straight out of the .html template rather
// than copied here, so these tests always exercise the code the server actually
// serves. A copy would drift and quietly stop catching regressions.

const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

const TEMPLATE_DIR = path.join(__dirname, "..", "..", "templates");

// the control logic lives in the one inline <script> block that defines
// controlButtons. the other blocks are figure/websocket plumbing.
function scriptBlock(template) {
    const src = fs.readFileSync(path.join(TEMPLATE_DIR, template), "utf8");
    const blocks = [...src.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((m) => m[1]);
    const found = blocks.filter((b) => b.includes("controlButtons"));
    if (found.length !== 1) {
        throw new Error(`${template}: expected one control script block, found ${found.length}`);
    }
    return found[0];
}

// every element id the template's control script reaches for
function referencedIds(template) {
    const matches = scriptBlock(template).matchAll(/getElementById\(['"]([^'"]+)['"]\)/g);
    return new Set([...matches].map((m) => m[1]));
}

const button = (id) => `<button type="button" id="${id}"></button>`;
const offButton = (id) => `<button type="button" id="${id}" disabled="disabled"></button>`;
const check = (id) => `<input type="checkbox" id="${id}">`;

// Minimal stand-ins for the control markup. Kept by hand rather than generated
// so that a template that grows a new control trips the fixture guard in
// controls.test.js instead of silently going untested.
const FIXTURES = {
    // the MMIRS page renders no mode selector, so the analyze url is built
    // differently there. mode: false reproduces that layout.
    "wfs.html": ({ mode = true } = {}) => `
        ${["focuscorrect", "comacorrect", "m1correct", "recenter"].map(offButton).join("")}
        ${["clear", "clearm1", "clearm2", "setgains", "analyze", "latest", "continuous"].map(button).join("")}
        <input id="m1gain" value="0.5">
        <input id="m2gain" value="1.0">
        ${["turbo", "connect", "spher"].map(check).join("")}
        <input id="datafile" value="">
        ${mode ? '<select id="mode"><option value="blue" selected></option></select>' : ""}
        <span id="datadir">/data/</span>
        <span id="zernikes"></span>`,
    "cwfs.html": () => `
        ${["focuscorrect", "comacorrect", "m1correct", "recenter"].map(offButton).join("")}
        ${["clear", "clearm1", "clearm2", "setgains", "analyze", "latest"].map(button).join("")}
        <input id="m1gain" value="0.5">
        <input id="m2gain" value="1.0">
        ${["connect", "spher"].map(check).join("")}
        <input id="datafile1" value="">
        <input id="datafile2" value="">
        <span id="datadir">/data/</span>
        <span id="zernikes"></span>`,
};

// let queued promise callbacks run. several rounds, because the handlers chain
// fetch promises a few deep (poll -> analyze -> corrections).
async function settle(rounds = 5) {
    for (let i = 0; i < rounds; i++) {
        await new Promise((resolve) => setImmediate(resolve));
    }
}

function mount(template, opts = {}) {
    const dom = new JSDOM(`<body>${FIXTURES[template](opts)}</body>`, { runScripts: "outside-only" });
    const win = dom.window;
    const doc = win.document;

    const calls = [];
    let files = opts.files || [];

    win.alert = () => {};
    win.fetch = (url, init) => {
        const u = String(url);
        calls.push({ url: u, path: u.split("?")[0], method: (init && init.method) || "GET" });
        const body = u.startsWith("files") ? JSON.stringify(files) : JSON.stringify("Z04 = 100 nm");
        return Promise.resolve({
            text: () => Promise.resolve(body),
            json: () => Promise.resolve(JSON.parse(body)),
        });
    };

    // virtual clock, so the 1s poll loop and the 5s running-flag reset fire in
    // the right order without the tests actually waiting on them
    let now = 0;
    let seq = 0;
    const timers = new Map();
    win.setTimeout = (fn, ms = 0) => {
        timers.set(++seq, { fn, due: now + ms, order: seq });
        return seq;
    };
    win.clearTimeout = (id) => timers.delete(id);

    async function advance(ms) {
        const target = now + ms;
        for (;;) {
            const due = [...timers.entries()]
                .filter(([, t]) => t.due <= target)
                .sort((a, b) => a[1].due - b[1].due || a[1].order - b[1].order);
            if (due.length === 0) break;
            const [id, timer] = due[0];
            timers.delete(id);
            now = timer.due;
            timer.fn();
            await settle();
        }
        now = target;
        await settle();
    }

    win.eval(scriptBlock(template));

    return {
        win,
        doc,
        calls,
        advance,
        settle,
        el: (id) => doc.getElementById(id),
        setFiles: (f) => {
            files = f;
        },
        urls: (p) => calls.filter((c) => c.path === p).map((c) => c.url),
        count: (p) => calls.filter((c) => c.path === p).length,
    };
}

module.exports = { mount, scriptBlock, referencedIds, FIXTURES };
