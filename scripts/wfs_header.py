#!/usr/bin/env python

import sys
from astropy.io import fits
from camsrv.header import update_header


if __name__ = "__main__":
    filename = sys.argv[1]
    hdu = fits.open(filename)
    newhdu = update_header(hdu)
    newhdu.writeto(filename, overwrite=True)
