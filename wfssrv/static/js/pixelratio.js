"use strict";

// mpl.js reads window.devicePixelRatio once, when a figure is made, and tells the server only when its websocket
// opens. Moving the window between a regular and a retina display changes the ratio, and the browser resizes the
// canvas to match, but the server keeps rendering at the old ratio, so the plots come out half (or double) size.
// This watches for the change and repeats what a fresh connection does: tell the server the new ratio and resize.
function watchPixelRatio(figures, win) {
    win = win || window;
    var current = win.devicePixelRatio || 1;

    function update() {
        var ratio = win.devicePixelRatio || 1;
        if (ratio === current) {
            return;
        }
        current = ratio;
        figures.forEach(function (fig) {
            if (fig.ws.readyState !== 1) {
                return;
            }
            fig.ratio = ratio;
            fig.send_message("set_device_pixel_ratio", { device_pixel_ratio: ratio });
            // the size mpl.js last gave the canvas, in css pixels. zero for a figure in a hidden tab, which gets
            // its size when the tab is shown.
            var width = parseFloat(fig.canvas.style.width) || 0;
            var height = parseFloat(fig.canvas.style.height) || 0;
            if (width > 0 && height > 0) {
                fig.canvas.setAttribute("width", Math.round(width * ratio));
                fig.canvas.setAttribute("height", Math.round(height * ratio));
                fig.request_resize(width, height);
            }
        });
    }

    // a resolution media query only matches the ratio it was made for, so make a new one after each change
    function listen() {
        var mq = win.matchMedia("(resolution: " + (win.devicePixelRatio || 1) + "dppx)");
        mq.addEventListener("change", function onchange() {
            mq.removeEventListener("change", onchange);
            update();
            listen();
        });
    }

    listen();
    // not every browser fires the media query change when the window moves to another display, so check on resize
    // as well. update() ignores anything that isn't a change of ratio.
    win.addEventListener("resize", update);
}

if (typeof module !== "undefined") {
    module.exports = { watchPixelRatio: watchPixelRatio };
}
