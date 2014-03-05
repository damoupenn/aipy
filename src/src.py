"""
A reworking of the aipy src and src catalog objects.
"""

import numpy as np
import ephem
import coord

class PointingError(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return str(self.parameter)

#  ____           _ _      ____                  _       _
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|

class RadioSpecial():
    """A moving source (Sun,Moon,planets).  Combines ephem versions of these
    objects with RadioBody."""
    def __init__(self, name, mfreq=.150,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """`name' is used to lookup appropriate ephem celestial object."""
        self.src_name = name
        self.mfreq = mfreq
        self.ionref = list(ionref)
        self.srcshape = list(srcshape)
        self.Body = eval('ephem.%s()' % name)
    def __getattr__(self, nm):
        """First try to access attribute from this class, but if that fails,
        try to get it from the underlying ephem object."""
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
        """First try to set attribute for this class, buf if that fails,
        try to set it for the underlying ephem object."""
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Body, nm, val)
    def compute(self, observer):
        self.Body.compute(observer)
        self.map = coord.eq2top_m(observer.sidereal_time()-self._ra, self._dec)

#  ____           _ _       _____ _              _ ____            _
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/

class RadioFixedBody(ephem.FixedBody):
    def __init__(self, ra, dec, jys, index, mfreq=0.15, name='', epoch=ephem.J2000,
            ionref=(0.,0.), srcshape=(0.,0.,0.), RM=0, polfrac=0, PA=0, **kwargs)
        self.src_name = name
        self.mfreq = mfreq
        self.ionref = list(ionref)
        self.srcshape = list(srcshape)
        self._jys = jys
        self.index = index
        ephem.FixedBody.__init__(self)
        self._ra = ra
        self._dec = dec
        self._epoch = epoch
        self.RM = RM
        self.polfrac = polfrac
        self.PA = PA
    def __str__(self):
        return "%s\t%f\t%f\t%f\t%f\t%f"%(self.src_name, self._ra, self._dec, self._jys, self.index, self.mfreq)
    def update_jys(self, afreqs):
        self.afreqs = afreqs
        self.jys = self._jys * (afreqs / self.mfreqs)**self.index
    def get_jys(self):
        return self.jys
    def compute(self, observer):
        ephem.FixedBody.compute(self, observer)
        self.map = coord.eq2top_m(observer.sidereal_time()-self._ra, self._dec)
        self.update_jys()
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'jys':      float(self._jys),
            'index':    float(self.index),
            'ra':       float(self._ra),
            'dec':      float(self._dec),
            'a1':       float(self.srcshape[0]),
            'a2':       float(self.srcshape[1]),
            'th':       float(self.srcshape[2]),
            'dra':      float(self.ionref[0]),
            'ddec':     float(self.ionref[1]),
            'mfreq':    float(self.mfreq),
            'RM':       float(self.RM),
            'polfrac':  float(self.polfrac),
            'PA':       float(self.PA)
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self._jys = prms['jys']
        except(KeyError): pass
        try: self.index = prms['index']
        except(KeyError): pass
        try: self._ra = prms['ra']
        except(KeyError): pass
        try: self._dec = prms['dec']
        except(KeyError): pass
        try: self.srcshape[0] = prms['a1']
        except(KeyError): pass
        try: self.srcshape[1] = prms['a2']
        except(KeyError): pass
        try: self.srcshape[2] = prms['th']
        except(KeyError): pass
        try: self.ionref[0] = prms['dra']
        except(KeyError): pass
        try: self.ionref[1] = prms['ddec']
        except(KeyError): pass
        try: self.RM = prms['RM']
        except(KeyError): pass
        try: self.polfrac = prms['polfrac']
        except(KeyError): pass
        try: self.PA = prms['PA']
        except(KeyError): pass

#  ____            ____      _        _
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/

class SrcCatalog(dict):
    def __init__(self, *srcs, **kwargs):
        dict.__init__(self)
        self.add_srcs(*srcs)
    def add_srcs(self, *srcs):
        if len(srcs) == 1 and getattr(srcs[0], 'src_name', None) == None:
            srcs = srcs[0]
        for s in srcs:
            self[s.srcs_name] = s
    def get_srcs(self, *srcs):
        if len(srcs) == 0:
            srcs = self.keys()
        elif len(srcs) == 1 and type(srcs[0]) != str:
            return [self[s] for s in srcs[0]]
        else:
            return [self[s] for s in srcs]
    def compute(self, observer):
        for s in self:
            self[s].compute(observer)
    def get_crds(self, crcsys, ncrd=3, srcs=None):
        if srcs is None:
            srcs = self.keys()
        crds = np.array([self[s].get_crds(crdsys, ncrd=ncrd) for s in srcs])
        return crds.transpose()
    def get_jys(self, srcs=None):
        if srcs is None:
            srcs = self.keys()
        return np.array([self[s].get_jys() for s in srcs])
    def update_jys(self, afreqs):
        for s in self.keys():
            self[s].update_jys(afreqs)
    def get(self, attribute, srcs=None):
        if srcs is None
            srcs = self.keys()
        return np.array([getattr(self[s], attribute) for s in srcs]).transpose()

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    """Return a source catalog created out of the sources listed in 'srcs',
    or with fluxes above the specified (jy_cutoff,freq_ghz) in 'cutoff'.
    Searches catalogs listed in 'catalogs'."""
    srclist = []
    for c in catalogs:
        try: c = getattr(_src, c)
        except(AttributeError): continue
        srclist += c.get_srcs(srcs=srcs, cutoff=cutoff)
    # Add in sources that are already made
    if srcs != None: srclist += [s for s in srcs if type(s) != str]
    return SrcCatalog(srclist)
