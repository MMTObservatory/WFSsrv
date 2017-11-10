#!/usr/bin/env python

import sys
import numpy as np
from astropy.io import fits
from camsrv.header import update_header


if __name__ == "__main__":
    filename = sys.argv[1]
    hdu = fits.open(filename)
    newhdu = update_header(hdu)

    # now check for identifying header information
    if isinstance(newhdu, fits.hdu.image.PrimaryHDU):
        header = newhdu.header
    elif isinstance(newhdu, list):
        header = None
        h = newhdu[1]
        header = h.header
        if header is None:
            raise ValueError("No PrimaryHDU found in HDU list.")
    else:
        raise ValueError("Must provide a PrimaryHDU object or an HDU list that contains one.")

    inst = None

    # check for MMIRS
    if 'WFSNAME' in header:
        if 'mmirs' in header['WFSNAME']:
            inst = "mmirs"

    # check for binospec
    if 'ORIGIN' in header:
        if 'Binospec' in header['ORIGIN']:
            inst = "binospec"
            # newhdu[1].data = np.flipud(newhdu[1].data)

    # check for F/9
    if 'OBSERVER' in header:
        if 'F/9 WFS' in header['OBSERVER']:
            inst = "f9"

    # check for F/5 (hecto)
    if inst is None and 'SEC' in header:  # mmirs has SEC in header as well which is why we check if inst isn't set
        if 'F5' in header['SEC']:
            inst = "f5"

    newhdu.writeto(filename, overwrite=True, output_verify="silentfix")

    print(inst)
