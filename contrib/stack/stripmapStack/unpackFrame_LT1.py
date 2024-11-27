#!/usr/bin/env python3

import isce
from isceobj.Sensor import createSensor
import shelve
import argparse
import glob
from isceobj.Util import Poly1D
from isceobj.Planet.AstronomicalHandbook import Const
import os


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='Unpack LT1 SLC data and store metadata in pickle file.')
    parser.add_argument('-i','--input', dest='h5dir', type=str,
            required=True, help='Input LT1 directory')
    parser.add_argument('-o', '--output', dest='slcdir', type=str,
            required=True, help='Output SLC directory')
    parser.add_argument('-d', '--deskew', dest='deskew', action='store_true',
            default=False, help='To read in for deskewing data later')
    parser.add_argument('-p', '--polarization', dest='polarization', type=str,
            default='HH', help='polarization in case if quad or full pol data exists. Deafult: HH')
    parser.add_argument('-orbit', '--orbitfile', dest='orbitfile', type=str, default=None, required=False,
                        help='Optional: directory with the precise orbit file for LT1 SLC. Default: None')
    return parser.parse_args()


def unpack(hdf5, slcname, deskew=False, polarization='HH', orbitfile=None):
    '''
    Unpack HDF5 to binary SLC file.
    '''

    tiffname = glob.glob(os.path.join(hdf5, '*.tiff'))[0]
    xmlname = glob.glob(os.path.join(hdf5, '*.meta.xml'))[0]

    
    if not os.path.isdir(slcname):
        os.mkdir(slcname)

    date = os.path.basename(slcname)
    obj = createSensor('LUTAN1')
    obj.configure()
    obj.xml = xmlname
    obj.tiff = tiffname

    if orbitfile:
        try:
            orbname = glob.glob(os.path.join(hdf5, '*ABSORBIT_SCIE.xml'), recursive=True)[0]
        except IndexError:
            try:
                orbname = glob.glob(os.path.join(hdf5, '*.txt'), recursive=True)[0]
            except IndexError:
                orbname = None
        obj.orbitFile = orbname
        print(obj.orbitFile)

    if deskew:
        obj.output = os.path.join(slcname, date+'_orig.slc')
    else:
        obj.output = os.path.join(slcname, date + '.slc')

    print(obj.xml)
    print(obj.tiff)
    print(obj.output)
    obj.extractImage()
    obj.frame.getImage().renderHdr()


    coeffs = obj.doppler_coeff
    dr = obj.frame.getInstrument().getRangePixelSize()
    r0 = obj.frame.getStartingRange()


    poly = Poly1D.Poly1D()
    poly.initPoly(order=1)
    poly.setCoeffs([0.0, 0.0])


    fpoly = Poly1D.Poly1D()
    fpoly.initPoly(order=1)
    fpoly.setCoeffs([0.0, 0.0])

    if deskew:
        pickName = os.path.join(slcname, 'original')
    else:
        pickName = os.path.join(slcname, 'data')
    with shelve.open(pickName) as db:
        db['frame'] = obj.frame
        db['doppler'] = poly
        db['fmrate'] = fpoly

    print("Doppler coefficients: ", poly._coeffs)
    print("FM rate coefficients: ", fpoly._coeffs)
    return obj

if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()
    if inps.slcdir.endswith('/'):
        inps.slcdir = inps.slcdir[:-1]

    if inps.h5dir.endswith('/'):
        inps.h5dir = inps.h5dir[:-1]

    obj = unpack(inps.h5dir, inps.slcdir,
                 deskew=inps.deskew,
                 polarization=inps.polarization,
                 orbitfile=inps.orbitfile)

