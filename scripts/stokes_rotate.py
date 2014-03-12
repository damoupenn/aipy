#!/usr/bin/env python
"""
Rotate linearly polarized data into Stokes' I,Q,U,V
"""
import aipy as a
import numpy as np
import optparse,sys,os

o = optparse.OptionParser()
o.set_usage('stokes_rotate.py *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


for uvfile in args:

    infile = uvfile
    outfile = infile+'S'

    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue

    uv = a.pol.UV(uvfile)
    DD = {}
    for (uvw,t,bl),d,f in uv.all(raw=True):
        plzn = uv.read_pol()
        if not bl in DD.keys(): DD[bl] = {}
        if not t in DD[bl].keys(): DD[bl][t] = {}
        if not plzn in DD[bl][t].keys():
            DD[bl][t][plzn] = np.ma.array(d,mask=f)
    del(uv)

    def mfunc(uv,p,d,f):
        uvw,t,bl = p
        plzn = uvi.read_pol()
        if plzn == 'xx':
            uvo.write_pol('I')
            return p,DD[bl][t]['xx'] + DD[bl][t]['yy'],f
        if plzn == 'xy':
            uvo.write_pol('Q')
            return p,DD[bl][t]['xx'] - DD[bl][t]['yy'],f
        if plzn == 'yx':
            uvo.write_pol('U')
            return p,DD[bl][t]['xy'] + DD[bl][t]['yx'],f
        if plzn == 'yy':
            uvo.write_pol('V')
            return p,-1.j*(DD[bl][t]['xy'] - DD[bl][t]['yx']),f

    uvi = a.pol.UV(infile)
    uvo = a.pol.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='XY --> STOKES \n')
