#!/usr/bin/env python

import sys
from astropy.io import fits
from camsrv.header import update_header


if __name__ == "__main__":
    filename = sys.argv[1]
    hdu = fits.open(filename)
    newhdu = update_header(hdu)
    newhdu.writeto(filename, overwrite=True, output_verify="silentfix")

    # now check for identifying header information
    if isinstance(newhdu, fits.hdu.image.PrimaryHDU):
        header = newhdu.header
    elif isinstance(newhdu, list):
        header = None
        for h in newhdu:
            if isinstance(h, fits.hdu.image.PrimaryHDU):
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

    # check for F/9
    if 'OBSERVER' in header:
        if 'F/9 WFS' in header['OBSERVER']:
            inst = "f9"

    # check for F/5 (hecto)
    if inst is None and 'SEC' in header:
        if 'F5' in header['SEC']:
            inst = "f5"

    print(inst)
