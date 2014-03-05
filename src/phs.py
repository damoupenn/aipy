"""
Module for representing antenna array geometry and for generating
phasing information.
"""
import ephem, math, numpy as n, coord, const, _cephes
from miriad import ij2bl, bl2ij

class PointingError(Exception):
    """An error to throw if a source is below the horizon."""
    def __init__(self, value): self.parameter = value
    def __str__(self): return str(self.parameter)

#  _   _ _   _ _ _ _           _____                 _   _
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/

def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Dec. 31, 1899."""
    return ephem.date(num - 2415020.)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Dec. 31, 1899) to Julian date."""
    return float(num + 2415020.)

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam:
    """Template for representing antenna beam pattern.  Beams also hold
    info about which frequencies are active (i.e. Antennas and AntennaArrays
    access frequencies through Beam)."""
    def __init__(self, freqs, **kwargs):
        """freqs = frequencies (in GHz) at bin centers across spectrum."""
        self.freqs = freqs
        self.chans = n.arange(self.freqs.size)
        self._update_afreqs()
    def _update_afreqs(self):
        self.afreqs = self.freqs.take(self.chans)
    def update(self):
        self._update_afreqs()
    def select_chans(self, active_chans=None):
        """Select only enumerated channels to use for future calculations."""
        if active_chans is None: active_chans = n.arange(self.freqs.size)
        self.chans = active_chans
        self.update()

#     _          _
#    / \   _ __ | |_ ___ _ __  _ __   __ _
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna:
    """Representation of physical attributes of individual antenna."""
    def __init__(self, x, y, z, beam, phsoff=[0.,0.], **kwargs):
        """x,y,z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object
        phsoff = polynomial phase vs. frequency.  Phs term that is linear
                 with freq is often called 'delay'."""
        self.pos = n.array((x,y,z), n.float64) # must be float64 for mir
        self.beam = beam
        self.__phsoff = phsoff
        self._update_phsoff()
        self.active_pol = None
    def select_chans(self, active_chans=None):
        """Select only the specified channels for use in future calculations."""
        self.beam.select_chans(active_chans)
        self.update()
    def set_active_pol(self, pol):
        assert(pol in 'xyIQUV') # We only do linear polarizations for now
        self.active_pol = pol
    def get_active_pol(self):
        if self.active_pol is None: raise RuntimeError('No active polarization set (use Antenna.set_active_pol)')
        return self.active_pol
    def _update_phsoff(self):
        self._phsoff = n.polyval(self.__phsoff, self.beam.afreqs)
    def update(self):
        self._update_phsoff()
    def phsoff(self):
        return self._phsoff
    def __iter__(self): return self.pos.__iter__()
    def __add__(self, a): return self.pos + a.pos
    __radd__ = __add__
    def __neg__(self): return -self.pos
    def __sub__(self, a): return self.pos - a.pos
    def __rsub__(self, a): return a.pos - self.pos

#     _                         _                    _   _
#    / \   _ __ _ __ __ _ _   _| |    ___   ___ __ _| |_(_) ___  _ __
#   / _ \ | '__| '__/ _` | | | | |   / _ \ / __/ _` | __| |/ _ \| '_ \
#  / ___ \| |  | | | (_| | |_| | |__| (_) | (_| (_| | |_| | (_) | | | |
# /_/   \_\_|  |_|  \__,_|\__, |_____\___/ \___\__,_|\__|_|\___/|_| |_|
#                         |___/

class ArrayLocation(ephem.Observer):
    """The location and time of an observation."""
    def __init__(self, location):
        """location = (lat,long,[elev]) of array"""
        ephem.Observer.__init__(self)
        self.pressure = 0
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
        self._update_eq2zen()
    def _update_eq2zen(self):
        self._eq2zen = coord.eq2top_m(0., self.lat)
    def update(self):
        self._update_eq2zen()
    def get_jultime(self):
        """Get current time as a Julian date."""
        return ephem2juldate(self.date)
    def set_jultime(self, t=None):
        """Set current time as a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set current time as derived from the ephem package.  Recalculates
        matrix for projecting baselines into current positions."""
        if t is None: t = ephem.now()
        self.date, self.epoch = t, t
        self._eq2now = coord.rot_m(-self.sidereal_time(), n.array([0.,0.,1.]))

#     _          _                            _
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/

class AntennaArray(ArrayLocation):
    """A collection of antennas, their spacings, and location/time of
    observations."""
    def __init__(self, location, ants, **kwargs):
        """ location = (lat,long,[elev]) of array
        ants = list of Antenna objects."""
        ArrayLocation.__init__(self, location=location)
        self.ants = ants
        self.active_pol = None
    def __iter__(self): return self.ants.__iter__()
    def __getitem__(self, *args): return self.ants.__getitem__(*args)
    def __setitem__(self, *args): return self.ants.__setitem__(*args)
    def __len__(self): return self.ants.__len__()
    def set_active_pol(self,pol):
        assert(pol in ('xx','xy','yx','yy','I','Q','U','V'))
        self.active_pol = pol
    def get_active_pol(self, split=False):
        if self.active_pol is None: raise RuntimeError('No active polarization set (use AntennaArray.set_active_pol)')
        if split: return self.active_pol[0], self.active_pol[-1]
        else: return self.active_pol
    def update(self):
        ArrayLocation.update(self)
        for a in self: a.update()
    def select_chans(self, active_chans=None):
        """Select which channels are used in computations.  Default is all."""
        for a in self: a.select_chans(active_chans)
        self.update()
    def ij2bl(self, i, j):
        """Convert baseline i,j (0 indexed) to Miriad's (i+1) << 8 | (j+1)
        indexing scheme."""
        return ij2bl(i,j)
    def bl2ij(self, bl):
        """Convert Miriad's (i+1) << 8 | (j+1) baseline indexing scheme to
        i,j (0 indexed)"""
        return bl2ij(bl)
    def bl_indices(self, auto=True, cross=True):
        """Return bl indices for baselines in the array."""
        if auto:
            if cross: return [self.ij2bl(i,j)
                for i in range(len(self))
                for j in range(i,len(self))]
            else: return [self.ij2bl(i,i) for i in range(len(self))]
        else:
            if cross: return [self.ij2bl(i,j)
                for i in range(len(self))
                for j in range(i+1,len(self))]
            else: return []
    def get_afreqs(self):
        """Return array of frequencies that are active for simulation."""
        return self[0].beam.afreqs
    def get_baseline(self, i, j, src='z'):
        """Return the baseline corresponding to i,j in various coordinate
        projections: src='e' for current equatorial, 'z' for zenith
        topocentric, 'r' for unrotated equatorial, or a RadioBody for
        projection toward that source."""
        bl = self[j] - self[i]
        if type(src) == str:
            if src == 'e': return n.dot(self._eq2now, bl)
            elif src == 'z': return n.dot(self._eq2zen, bl)
            elif src == 'r': return bl
            else: raise ValueError('Unrecognized source:' + src)
        try:
            if src.alt < 0:
                raise PointingError('%s below horizon' % src.src_name)
            m = src.map
        except(AttributeError):
            ra,dec = coord.eq2radec(src)
            m = coord.eq2top_m(self.sidereal_time() - ra, dec)
        return n.dot(m, bl).transpose()
    def get_phs_offset(self, i, j):
        """Return the frequency-dependent phase offset of baseline i,j."""
        pi,pj = self.get_active_pol(split=True)
        self[i].set_active_pol(pi)
        self[j].set_active_pol(pj)
        self.set_active_pol(pi+pj)
        return self[j].phsoff() - self[i].phsoff()
    def gen_uvw(self, i, j, src='z', w_only=False):
        """Compute uvw coordinates of baseline relative to provided RadioBody,
        or 'z' for zenith uvw coordinates.  If w_only is True, only w (instead
        of (u,v,w) will be returned)."""
        x,y,z = self.get_baseline(i,j, src=src)
        afreqs = self.get_afreqs()
        afreqs = n.reshape(afreqs, (1,afreqs.size))
        if len(x.shape) == 0:
            if w_only: return z*afreqs
            else: return n.array([x*afreqs, y*afreqs, z*afreqs])
        #afreqs = n.reshape(afreqs, (1,afreqs.size))
        x.shape += (1,); y.shape += (1,); z.shape += (1,)
        if w_only: return n.dot(z,afreqs)
        else: return n.array([n.dot(x,afreqs), n.dot(y,afreqs), n.dot(z,afreqs)])
    def gen_phs(self, src, i, j, mfreq=.150, ionref=None, srcshape=None,
            resolve_src=False, nocal=False):
        """Return phasing that is multiplied to data to point to src."""
        if ionref is None:
            try: ionref = src.ionref
            except(AttributeError): pass
        if not ionref is None or resolve_src: u,v,w = self.gen_uvw(i,j,src=src)
        else: w = self.gen_uvw(i,j,src=src, w_only=True)
        if not ionref is None: w += self.refract(u, v, mfreq=mfreq, ionref=ionref)
        phs = n.exp(-1j*2*n.pi*(w))
        if not nocal: phs *= n.exp(-1j*2*n.pi*(self.get_phs_offset(i,j)))
        if resolve_src:
            if srcshape is None:
                try: srcshape = src.srcshape
                except(AttributeError): pass
            if not srcshape is None: phs *= self.resolve_src(u, v, srcshape=srcshape)
        return phs.squeeze()
    def resolve_src(self, u, v, srcshape=(0,0,0)):
        """Adjust amplitudes to reflect resolution effects for a uniform
        elliptical disk characterized by srcshape:
        srcshape = (a1,a2,th) where a1,a2 are angular sizes along the
            semimajor, semiminor axes, and th is the angle (in radians) of
            the semimajor axis from E."""
        a1,a2,th = srcshape
        try:
            if len(u.shape) > len(a1.shape):
                a1.shape += (1,); a2.shape += (1,); th.shape += (1,)
        except(AttributeError): pass
        ru = a1 * (u*n.cos(th) - v*n.sin(th))
        rv = a2 * (u*n.sin(th) + v*n.cos(th))
        x = 2 * n.pi * n.sqrt(ru**2 + rv**2)
        # Use first Bessel function of the first kind (J_1)
        return n.where(x == 0, 1, 2 * _cephes.j1(x)/x).squeeze()
    def refract(self, u_sf, v_sf, mfreq=.150, ionref=(0.,0.)):
        """Calibrate a frequency-dependent source offset by scaling measured
        offsets at a given frequency.  Generates dw, a change in the
        projection of a baseline towards that source, which can be used to
        fix the computed phase of that source.
        ionref = (dra, ddec) where dra, ddec are angle offsets (in radians)
            of sources along ra/dec axes at the specified mfreq.
        u_sf,v_sf = u,v components of baseline, used to compute the
            change in w given angle offsets and the small angle approx.  Should
            be numpy arrays with sources (s) along the 1st axis and
            freqs (f) along the 2nd."""
        dra,ddec = ionref
        s,f = u_sf.shape
        try: dra.shape = (s,1)
        except(AttributeError): pass
        try: ddec.shape = (s,1)
        except(AttributeError): pass
        try: mfreq.shape = (s,1)
        except(AttributeError): pass
        f2 = self.get_afreqs()**2 ; f2.shape = (1, f)
        return (dra*u_sf + ddec*v_sf) * mfreq**2 / f2
    def phs2src(self, data, src, i, j, mfreq=.150, ionref=None, srcshape=None, nocal=False):
        """Apply phasing to zenith-phased data to point to src."""
        return data * self.gen_phs(src, i, j,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=None, nocal=nocal)
    def unphs2src(self,data,src, i, j, mfreq=.150, ionref=None, srcshape=None, nocal=False):
        """Remove phasing from src-phased data to point to zenith."""
        return data / self.gen_phs(src, i, j,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=None, nocal=nocal)
