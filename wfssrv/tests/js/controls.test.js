"use strict";

// Tests for the control logic in the wfssrv templates, in particular the
// continuous WFS loop. See harness.js: the javascript is read out of the
// templates, so these tests track what the server actually serves.

const test = require("node:test");
const assert = require("node:assert");
const { mount, referencedIds } = require("./harness.js");

// a new frame has landed on top of the one currently loaded
const NEW_FRAME = ["wfs_0002.fits", "wfs_0001.fits"];

function startContinuous(page, { connect = false, turbo = false } = {}) {
    page.el("connect").checked = connect;
    page.el("turbo").checked = turbo;
    page.el("datafile").value = "wfs_0001.fits";
    page.el("continuous").click();
}

test("continuous mode analyzes a new frame even though the controls are disabled", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    await page.advance(1500);

    // the loop disables the controls to keep the operator's hands off them, so
    // it must not rely on the buttons being clickable. this is the regression
    // that broke continuous mode when jquery was removed.
    assert.equal(page.el("analyze").disabled, true, "analyze should be disabled during continuous");
    assert.equal(page.count("analyze"), 1, "a new frame should have been analyzed");
    assert.match(page.urls("analyze")[0], /fitsfile=\/data\/wfs_0002\.fits/);
});

test("continuous mode applies focus, coma and m1 corrections when connected", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page, { connect: true });
    await page.advance(1500);

    assert.equal(page.count("focuscorrect"), 1);
    assert.equal(page.count("comacorrect"), 1);
    assert.equal(page.count("m1correct"), 1);
});

test("continuous mode sends no corrections when not connected", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page, { connect: false });
    await page.advance(1500);

    assert.equal(page.count("analyze"), 1, "it should still analyze");
    assert.equal(page.count("focuscorrect") + page.count("comacorrect") + page.count("m1correct"), 0);
});

test("continuous mode keeps polling but does not re-analyze the same frame", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    await page.advance(20000);

    assert.ok(page.count("files") > 2, "the loop should keep polling for new frames");
    assert.equal(page.count("analyze"), 1, "the same frame should only be analyzed once");
});

test("continuous mode picks up each new frame as it arrives", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    await page.advance(10000);
    assert.equal(page.count("analyze"), 1);

    page.setFiles(["wfs_0003.fits", ...NEW_FRAME]);
    await page.advance(10000);

    assert.equal(page.count("analyze"), 2);
    assert.match(page.urls("analyze")[1], /fitsfile=\/data\/wfs_0003\.fits/);
});

test("the MMIRS layout has no mode selector and omits mode from the analyze url", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME, mode: false });
    assert.equal(page.el("mode"), null, "MMIRS renders no mode selector");
    startContinuous(page);
    await page.advance(1500);

    const url = page.urls("analyze")[0];
    assert.doesNotMatch(url, /mode=/);
    assert.match(url, /fitsfile=\/data\/wfs_0002\.fits/);
});

test("a layout with a mode selector passes the mode through", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    await page.advance(1500);

    assert.match(page.urls("analyze")[0], /mode=blue/);
});

test("starting continuous mode backs the gains off, unless turbo is set", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    assert.equal(Number(page.el("m1gain").value), 0.2);
    assert.equal(Number(page.el("m2gain").value), 0.5);

    const turbo = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(turbo, { turbo: true });
    assert.equal(Number(turbo.el("m1gain").value), 0.5);
    assert.equal(Number(turbo.el("m2gain").value), 1.0);
});

test("stopping continuous mode restores the gains and re-enables the controls", async () => {
    const page = mount("wfs.html", { files: NEW_FRAME });
    startContinuous(page);
    await page.advance(1500);

    page.el("continuous").click();
    await page.settle();

    assert.equal(Number(page.el("m1gain").value), 0.5);
    assert.equal(Number(page.el("m2gain").value), 1.0);
    assert.equal(page.el("analyze").disabled, false);
    assert.equal(page.count("m1gain"), 2, "gains should be pushed on start and on stop");

    const analyses = page.count("analyze");
    await page.advance(20000);
    assert.equal(page.count("analyze"), analyses, "the loop should be stopped");
});

test("the cwfs latest button loads the two newest frames and analyzes them", async () => {
    const page = mount("cwfs.html", { files: ["cwfs_0002.fits", "cwfs_0001.fits"] });
    page.el("latest").click();
    await page.settle();

    assert.equal(page.el("datafile1").value, "cwfs_0002.fits");
    assert.equal(page.el("datafile2").value, "cwfs_0001.fits");
    assert.equal(page.count("analyze"), 1);
});

// The fixtures in harness.js are maintained by hand. This keeps them honest: if
// a template grows or renames a control, this fails and tells you to update the
// fixture rather than quietly testing less than it appears to.
for (const template of ["wfs.html", "cwfs.html"]) {
    test(`the ${template} fixture provides every element the template uses`, () => {
        const page = mount(template);
        for (const id of referencedIds(template)) {
            assert.ok(page.el(id), `${template}: fixture is missing #${id}`);
        }
    });
}
