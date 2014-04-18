#! /usr/bin/env python

"""
    Rotate linearly polarized data into Stokes' IQUV
"""

import aipy as a
import numpy as np
import optparse, sys, os

o = optparse.OptionParser()
o.set_usage('StokesRotate.py *.uv')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])

for uvfile in args:
    infile = uvfile
    outfile = infile+'S'

    print infile, '-->', outfile

    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue

    uv = a.pol.UV(uvfile)
    DD = {}
    for (uvw,t,bl),d,f in uv.all(raw=True):
        pol = uv.read_pol()
        if not bl in DD.keys(): DD[bl] = {}
        if not t in DD[bl].keys(): DD[bl][t] = {}
        if not pol in DD[bl][t].keys():
            DD[bl][t][pol] = np.ma.array(d, mask=f)
    del uv

    uvi = a.pol.UV(infile)
    uvo = a.pol.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    #this is just being more explicit about header info than uv.all()
    while True:
        try:
            p,d,f = uvi.read(raw=True)
            uvw,t,bl = p
            pol = uvi.read_pol()
            try:
                if pol == 'xx':
                    uvo.write_pol('I')
                    f = DD[bl][t]['xx'].mask | DD[bl][t]['yy'].mask
                    d = DD[bl][t]['xx'] + DD[bl][t]['yy']
                elif pol == 'xy':
                    uvo.write_pol('Q')
                    f = DD[bl][t]['xx'].mask | DD[bl][t]['yy'].mask
                    d = DD[bl][t]['xx'] - DD[bl][t]['yy']
                elif pol == 'yx':
                    uvo.write_pol('U')
                    f = DD[bl][t]['xy'].mask | DD[bl][t]['yx'].mask
                    d = DD[bl][t]['xy'] + DD[bl][t]['yx']
                elif pol == 'yy':
                    uvo.write_pol('V')
                    f = DD[bl][t]['xy'].mask | DD[bl][t]['yx'].mask
                    d = -1.j*(DD[bl][t]['xy'] - DD[bl][t]['yx'])
                uvo.write(p,d,f)
            except(KeyError):
                pass
        except(IOError):
            break

    uvo._wrhd('history', uvo['history'] + "XY --> STOKES \n")

    del uvi, uvo
