"""
sruth_vmadcp_core.py
Processing engine for UHDAS / Teledyne Ocean Surveyor (os75bb / os75nb) cruise folders.

Two data paths, both producing the same Section object:
  1. load_codas()   - reads CODAS output already in proc/<sonar>/contour/<sonar>.nc
  2. process_raw()  - single-ping processing from raw/os75/*.raw (+ .raw.log),
                      raw/seapath/*.sea and raw/gyro/*.hdg

Post-processing applied to a Section (either path):
  apply_calibration() - rotation (deg) + amplitude scale factor on ship-relative velocity
  apply_editing()     - PG threshold, CODAS flags, depth/time windows, speed & spike limits

Only numpy / scipy / matplotlib are required. netCDF4 is optional
(scipy.io.netcdf_file is used as fallback for reading/writing classic netCDF).
"""
import os, re, glob, struct, math, datetime as dt
import numpy as np

YEARBASE_DEFAULT = 2026
BAD_VEL = -32768


# ----------------------------------------------------------------------------- utils
def dday_to_datetime(dday, yearbase):
    base = dt.datetime(int(yearbase), 1, 1)
    return [base + dt.timedelta(days=float(d)) for d in np.atleast_1d(dday)]


def uv_to_spd_dir(u, v):
    spd = np.hypot(u, v)
    direc = (np.degrees(np.arctan2(u, v)) + 360.0) % 360.0   # direction current flows TO
    return spd, direc


def rotate_uv(u, v, ang_deg):
    """Rotate vectors clockwise by ang_deg (CODAS 'rotate' sign convention:
    positive angle turns the velocity clockwise, i.e. same as adding to heading)."""
    a = np.radians(ang_deg)
    ca, sa = np.cos(a), np.sin(a)
    return u * ca + v * sa, -u * sa + v * ca


def _read_text(path):
    try:
        with open(path, 'r', errors='replace') as f:
            return f.read()
    except Exception:
        return ''


# ----------------------------------------------------------------------------- folder discovery
class CruiseFolder:
    """Inspects a UHDAS cruise directory and reports what's available."""

    def __init__(self, root):
        self.root = os.path.abspath(root)
        if not os.path.isdir(self.root):
            raise FileNotFoundError(root)
        self.name = os.path.basename(self.root.rstrip('\\/'))
        self.proc_sonars = sorted(d for d in os.listdir(self._p('proc'))
                                  if os.path.isdir(self._p('proc', d))) if os.path.isdir(self._p('proc')) else []
        self.raw_adcp_dirs = []
        rawdir = self._p('raw')
        self.raw_subdirs = sorted(os.listdir(rawdir)) if os.path.isdir(rawdir) else []
        for d in self.raw_subdirs:
            if glob.glob(os.path.join(rawdir, d, '*.raw')):
                self.raw_adcp_dirs.append(d)
        self.config = self._read_proc_config()
        self.cals = parse_cals(_read_text(self._p('reports', 'ADCP_cals.txt')))

    def _p(self, *a):
        return os.path.join(self.root, *a)

    def nc_file(self, sonar):
        f = self._p('proc', sonar, 'contour', sonar + '.nc')
        return f if os.path.exists(f) else None

    def raw_files(self, inst='os75'):
        return sorted(glob.glob(self._p('raw', inst, '*.raw')))

    def nav_files(self, sub, ext):
        return sorted(glob.glob(self._p('raw', sub, '*.' + ext)))

    def _read_proc_config(self):
        cfg = {}
        cands = glob.glob(self._p('raw', 'config', '*_proc.toml'))
        if not cands:
            cands = glob.glob(self._p('proc', '*', 'config', '*_proc.toml'))
        if not cands:
            return cfg
        txt = _read_text(cands[0])
        cfg['_file'] = cands[0]
        for line in txt.splitlines():
            line = line.split('#')[0].strip()
            m = re.match(r'^(\w+)\s*=\s*(.+)$', line)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if v.startswith('{'):
                d = {}
                for kk, vv in re.findall(r'(\w+)\s*=\s*("?[^,}"]*"?)', v):
                    d[kk] = _cast(vv)
                cfg[k] = d
            else:
                cfg[k] = _cast(v)
        return cfg

    def cfg_get(self, key, inst='os75', default=None):
        v = self.config.get(key, default)
        if isinstance(v, dict):
            return v.get(inst, v.get(inst[:4], default))
        return v

    def summary(self):
        lines = [f'Cruise folder : {self.root}',
                 f'Cruise        : {self.config.get("cruiseid", self.name)}   Ship: {self.config.get("shipname", "?")}',
                 f'CODAS sonars  : {", ".join(self.proc_sonars) or "none"}',
                 f'Raw ADCP dirs : {", ".join(self.raw_adcp_dirs) or "none"}',
                 f'Raw subdirs   : {", ".join(self.raw_subdirs)}']
        for s in self.proc_sonars:
            lines.append(f'  {s}: netCDF {"found" if self.nc_file(s) else "MISSING"}')
        for inst in self.raw_adcp_dirs:
            lines.append(f'  raw/{inst}: {len(self.raw_files(inst))} .raw files')
        for sub, ext in (('seapath', 'sea'), ('gyro', 'hdg'), ('cnav', 'gps')):
            n = len(self.nav_files(sub, ext))
            if n:
                lines.append(f'  raw/{sub}: {n} .{ext} files')
        if self.config:
            lines.append('Processing config (' + os.path.basename(self.config.get('_file', '')) + '):')
            for k in ('h_align', 'ducer_depth', 'enslength', 'xducer_dx', 'xducer_dy', 'scalefactor',
                      'hdg_inst', 'hcorr_inst', 'pos_inst'):
                if k in self.config:
                    lines.append(f'  {k:12s} = {self.config[k]}')
        if self.cals:
            lines.append('Calibrations found in reports/ADCP_cals.txt:')
            for s, d in self.cals.items():
                lines.append('  ' + s + ': ' + ', '.join(f'{k}={v}' for k, v in d.items()))
        return '\n'.join(lines)


def _cast(v):
    v = v.strip().strip('"').strip("'")
    for t in (int, float):
        try:
            return t(v)
        except ValueError:
            pass
    return v


def parse_cals(txt):
    """Pull BT/WT median amplitude & phase and DXDY guesses from ADCP_cals.txt."""
    out = {}
    for line in txt.splitlines():
        m = re.match(r'cal:\s+(\S+)\s+(BT|WT)\s+(amplitude|phase)\s+([-\d.]+)\s+([-\d.]+)', line)
        if m:
            s, kind, q, med = m.group(1), m.group(2), m.group(3), float(m.group(4))
            out.setdefault(s, {})[f'{kind}_{q[:3]}'] = med
        m = re.match(r'cal:\s+(\S+)\s+DXDY\s+xducer_(dx|dy)\s*=\s*([-\d.]+)', line)
        if m:
            out.setdefault(m.group(1), {})[m.group(2)] = float(m.group(3))
    return out


# ----------------------------------------------------------------------------- Section container
class Section:
    """Gridded ensemble data: time x depth_cell. Velocities are earth-referenced (absolute)."""

    def __init__(self, **kw):
        self.name = kw.get('name', '')
        self.source = kw.get('source', '')
        self.yearbase = kw.get('yearbase', YEARBASE_DEFAULT)
        self.time = np.asarray(kw['time'], float)                  # decimal day
        self.lon = np.asarray(kw['lon'], float)
        self.lat = np.asarray(kw['lat'], float)
        self.depth = np.asarray(kw['depth'], float)                # (nt, nz)
        self.u = np.asarray(kw['u'], float)
        self.v = np.asarray(kw['v'], float)
        nt, nz = self.u.shape
        nanz = np.full((nt, nz), np.nan)
        self.w = np.asarray(kw.get('w', nanz), float)
        self.e = np.asarray(kw.get('e', nanz), float)              # error velocity
        self.amp = np.asarray(kw.get('amp', nanz), float)
        self.pg = np.asarray(kw.get('pg', np.full((nt, nz), 100.0)), float)
        self.pflag = np.asarray(kw.get('pflag', np.zeros((nt, nz))), float)
        self.uship = np.asarray(kw.get('uship', np.full(nt, np.nan)), float)
        self.vship = np.asarray(kw.get('vship', np.full(nt, np.nan)), float)
        self.heading = np.asarray(kw.get('heading', np.full(nt, np.nan)), float)
        self.tr_temp = np.asarray(kw.get('tr_temp', np.full(nt, np.nan)), float)
        self.num_pings = np.asarray(kw.get('num_pings', np.full(nt, np.nan)), float)
        self.bt_depth = np.asarray(kw.get('bt_depth', np.full(nt, np.nan)), float)
        self.meta = dict(kw.get('meta', {}))
        self.log = list(kw.get('log', []))
        self.mask = ~np.isfinite(self.u)          # True = rejected

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def mid_time(self):
        """Ensemble mid-times (CODAS stamps ensembles at their end)."""
        if self.meta.get('time_is_end'):
            es = self.meta.get('ensemble_seconds')
            es = float(np.nanmedian(es)) if es is not None else 300.0
            return self.time - 0.5 * es / 86400.0
        return self.time

    # derived
    @property
    def spd(self):
        return np.hypot(self.um, self.vm)

    @property
    def um(self):
        return np.where(self.mask, np.nan, self.u)

    @property
    def vm(self):
        return np.where(self.mask, np.nan, self.v)

    @property
    def direction(self):
        return uv_to_spd_dir(self.um, self.vm)[1]

    def distance_km(self):
        lat = np.radians(self.lat); lon = np.radians(self.lon)
        d = np.zeros_like(lat)
        ok = np.isfinite(lat) & np.isfinite(lon)
        idx = np.where(ok)[0]
        if len(idx) > 1:
            dl = np.diff(lat[idx]); dn = np.diff(lon[idx])
            a = np.sin(dl / 2) ** 2 + np.cos(lat[idx][:-1]) * np.cos(lat[idx][1:]) * np.sin(dn / 2) ** 2
            seg = 2 * 6371.0 * np.arcsin(np.sqrt(a))
            d[idx] = np.concatenate([[0], np.cumsum(seg)])
            d[~ok] = np.nan
        return d

    def layer_mean(self, zmin, zmax):
        sel = (self.depth >= zmin) & (self.depth <= zmax)
        u = np.where(sel, self.um, np.nan); v = np.where(sel, self.vm, np.nan)
        with np.errstate(all='ignore'):
            return np.nanmean(u, axis=1), np.nanmean(v, axis=1)

    def mean_profile(self):
        z = np.nanmedian(self.depth, axis=0)
        with np.errstate(all='ignore'):
            return z, np.nanmean(self.um, axis=0), np.nanmean(self.vm, axis=0), np.nanstd(self.um, axis=0), np.nanstd(self.vm, axis=0)

    def stats_text(self):
        good = (~self.mask).sum()
        tot = self.mask.size
        s = [f'{self.name}  [{self.source}]',
             f'  profiles: {len(self.time)}   bins: {self.u.shape[1]}   good cells: {good}/{tot} ({100 * good / max(tot, 1):.1f}%)',
             f'  time: {self.time[0]:.4f} - {self.time[-1]:.4f} dday   '
             f'({dday_to_datetime(self.time[0], self.yearbase)[0]:%Y-%m-%d %H:%M} to '
             f'{dday_to_datetime(self.time[-1], self.yearbase)[0]:%Y-%m-%d %H:%M} UTC)',
             f'  lon: {np.nanmin(self.lon):.4f} to {np.nanmax(self.lon):.4f}   lat: {np.nanmin(self.lat):.4f} to {np.nanmax(self.lat):.4f}',
             f'  depth: {np.nanmin(self.depth):.1f} - {np.nanmax(self.depth):.1f} m']
        with np.errstate(all='ignore'):
            sp = self.spd
            if np.isfinite(sp).any():
                s.append(f'  speed: mean {np.nanmean(sp):.3f}  max {np.nanmax(sp):.3f} m/s')
            deepest = np.where(~self.mask, self.depth, np.nan)
            if np.isfinite(deepest).any():
                s.append(f'  deepest good bin (median over profiles): {np.nanmedian(np.nanmax(deepest, axis=1)):.1f} m')
        return '\n'.join(s)


# ----------------------------------------------------------------------------- CODAS netCDF
def _open_nc(path):
    try:
        import netCDF4
        ds = netCDF4.Dataset(path)
        ds.set_auto_mask(True)

        def get(k):
            if k not in ds.variables:
                return None
            a = ds.variables[k][:]
            return np.ma.filled(np.ma.asarray(a).astype(float), np.nan)
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
        return get, attrs, ds.close
    except ImportError:
        from scipy.io import netcdf_file
        ds = netcdf_file(path, 'r', mmap=False, maskandscale=True)

        def get(k):
            if k not in ds.variables:
                return None
            a = np.array(ds.variables[k][:], dtype=float)
            fv = getattr(ds.variables[k], '_FillValue', None)
            if fv is not None:
                a[a == fv] = np.nan
            a[np.abs(a) > 1e30] = np.nan
            return a
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in ds._attributes.items()}
        return get, attrs, ds.close


def load_codas(cf, sonar):
    path = cf.nc_file(sonar)
    if path is None:
        raise FileNotFoundError(f'No CODAS netCDF for {sonar}')
    get, attrs, close = _open_nc(path)
    try:
        yb = int(attrs.get('yearbase', cf.config.get('yearbase', YEARBASE_DEFAULT)))
        sec = Section(name=sonar, source='CODAS netCDF', yearbase=yb,
                      time=get('time'), lon=get('lon'), lat=get('lat'), depth=get('depth'),
                      u=get('u'), v=get('v'), amp=get('amp'), pg=get('pg'),
                      pflag=np.nan_to_num(get('pflag')), uship=get('uship'), vship=get('vship'),
                      heading=get('heading'), tr_temp=get('tr_temp'), num_pings=get('num_pings'),
                      meta={'file': path, 'ensemble_seconds': get('ensemble_seconds'), 'time_is_end': True,
                            'transducer_depth': get('transducer_depth'),
                            'depth_bin_length': get('depth_bin_length'),
                            'transducer_orientation': get('transducer_orientation')})
    finally:
        close()
    # bottom depth from allbins_other.mat if present (LGB/MAB not always reliable; use for info)
    try:
        from scipy.io import loadmat
        m = loadmat(os.path.join(cf.root, 'proc', sonar, 'contour', 'allbins_other.mat'))
        lgb = m['LGB'].ravel().astype(float)
        if len(lgb) == len(sec.time):
            dep = sec.depth
            idx = np.clip(lgb.astype(int) - 1, 0, dep.shape[1] - 1)
            bt = dep[np.arange(len(idx)), idx]
            bt[~np.isfinite(lgb) | (lgb <= 0)] = np.nan
            sec.meta['last_good_bin_depth'] = bt
    except Exception:
        pass
    # u/v in the short netCDF are already masked where CODAS flagged them. Keep a copy of the
    # full-range (unflagged) velocities is not possible from the short form; pflag retained.
    sec.log.append(f'Loaded {os.path.basename(path)}: {len(sec.time)} profiles, {sec.u.shape[1]} bins')
    return sec


# ----------------------------------------------------------------------------- NMEA / nav readers
def _nmea_ok(s):
    return '*' in s


def read_seapath(files, progress=None):
    """Returns dict of arrays: t (dday PC clock), lat, lon (from INGGA/GPGGA);
    th, hdg (INHDT or PSXN,23); tp, pitch, roll (PSXN,23)."""
    gt, glat, glon = [], [], []
    ht, hh = [], []
    at, ap, ar = [], [], []
    for k, f in enumerate(files):
        if progress:
            progress(f'reading {os.path.basename(f)}', k / max(len(files), 1))
        tcur = None
        with open(f, 'r', errors='replace') as fh:
            for line in fh:
                if line.startswith('$PYRTM'):
                    p = line.split(',')
                    try:
                        tcur = float(p[2])
                    except (IndexError, ValueError):
                        tcur = None
                    continue
                if tcur is None:
                    continue
                if line[3:6] == 'GGA':
                    p = line.split(',')
                    try:
                        if p[6] in ('0', ''):
                            continue
                        la = float(p[2][:2]) + float(p[2][2:]) / 60.0
                        if p[3] == 'S':
                            la = -la
                        lo = float(p[4][:3]) + float(p[4][3:]) / 60.0
                        if p[5] == 'W':
                            lo = -lo
                        gt.append(tcur); glat.append(la); glon.append(lo)
                    except (IndexError, ValueError):
                        pass
                elif line[3:6] == 'HDT':
                    try:
                        hh.append(float(line.split(',')[1])); ht.append(tcur)
                    except (IndexError, ValueError):
                        pass
                elif line.startswith('$PSXN,23'):
                    p = line.split(',')
                    try:
                        ar.append(float(p[2])); ap.append(float(p[3])); at.append(tcur)
                    except (IndexError, ValueError):
                        pass
                tcur = None
    return dict(t=np.array(gt), lat=np.array(glat), lon=np.array(glon),
                th=np.array(ht), hdg=np.array(hh),
                tp=np.array(at), pitch=np.array(ap), roll=np.array(ar))


def read_gyro(files):
    t, h = [], []
    for f in files:
        tcur = None
        with open(f, 'r', errors='replace') as fh:
            for line in fh:
                if line.startswith('$PYRTM'):
                    try:
                        tcur = float(line.split(',')[2])
                    except (IndexError, ValueError):
                        tcur = None
                elif tcur is not None and line[3:6] == 'HDT':
                    try:
                        h.append(float(line.split(',')[1])); t.append(tcur)
                    except (IndexError, ValueError):
                        pass
                    tcur = None
    return dict(t=np.array(t), hdg=np.array(h))


def interp_angle(t, tp, ang):
    if len(tp) < 2:
        return np.full(len(t), np.nan)
    o = np.argsort(tp)
    tp, ang = tp[o], ang[o]
    r = np.unwrap(np.radians(ang))
    out = np.degrees(np.interp(t, tp, r, left=np.nan, right=np.nan)) % 360.0
    return out


def angdiff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


# ----------------------------------------------------------------------------- PD0 reader (Ocean Surveyor interleaved bb/nb)
# OS interleaved: one ping type under IDs 0x0000/0x0080/0x0100/0x0200/0x0300, the other under 0x0001/0x0081/...;
# the type is read from fixed-leader byte 16 (1 = bb, 10 = nb). BT = 0x0600.
PINGTYPE_IDS = {'bb': dict(fixed=0x0001, var=0x0081, vel=0x0101, cor=0x0201, amp=0x0301),
                'nb': dict(fixed=0x0000, var=0x0080, vel=0x0100, cor=0x0200, amp=0x0300)}


def _parse_fixed(b, o):
    return dict(ncell=b[o + 9], npings=struct.unpack_from('<H', b, o + 10)[0],
                cell=struct.unpack_from('<H', b, o + 12)[0] / 100.0,
                blank=struct.unpack_from('<H', b, o + 14)[0] / 100.0,
                coord=b[o + 25],
                ea=struct.unpack_from('<h', b, o + 26)[0] / 100.0,
                bin1=struct.unpack_from('<H', b, o + 32)[0] / 100.0,
                xmit=struct.unpack_from('<H', b, o + 34)[0] / 100.0,
                sysconf=struct.unpack_from('<H', b, o + 4)[0])


def read_pd0_file(path, pingtype, maxcell=128):
    """Return list-free dict of arrays for one .raw file and one ping type."""
    b = open(path, 'rb').read()
    logf = path + '.log'
    dday = None
    if os.path.exists(logf):
        try:
            L = np.loadtxt(logf, ndmin=2)
            dday = dict(zip(L[:, 1].astype(np.int64), L[:, 4]))   # byte offset -> PC dday
        except Exception:
            dday = None
    recs = []
    n = len(b)
    i = 0
    while i < n - 6:
        if b[i] != 0x7f or b[i + 1] != 0x7f:
            j = b.find(b'\x7f\x7f', i + 1)
            if j < 0:
                break
            i = j
            continue
        nbytes = struct.unpack_from('<H', b, i + 2)[0]
        if i + nbytes + 2 > n:
            break
        chk = sum(b[i:i + nbytes]) & 0xFFFF
        if chk != struct.unpack_from('<H', b, i + nbytes)[0]:
            i += 1
            continue
        nd = b[i + 5]
        offs = struct.unpack_from('<' + 'H' * nd, b, i + 6)
        blocks = {struct.unpack_from('<H', b, i + o)[0]: i + o for o in offs}
        # Identify the ping type from the fixed-leader mode byte (1 = broadband, 10 = narrowband);
        # with only one ping type enabled the OS writes it under the 0x0000 IDs regardless of type.
        ids = None
        for suf in (0x0000, 0x0001):
            if suf in blocks and (0x0100 | suf) in blocks:
                mode = b[blocks[suf] + 16]
                if ('bb' if mode == 1 else 'nb') == pingtype:
                    ids = dict(fixed=suf, var=0x0080 | suf, vel=0x0100 | suf, cor=0x0200 | suf, amp=0x0300 | suf)
                    break
        if ids is not None:
            fl = _parse_fixed(b, blocks[ids['fixed']])
            nc = fl['ncell']
            vo = blocks[ids['vel']] + 2
            vel = np.frombuffer(b, '<i2', nc * 4, vo).reshape(nc, 4).astype(float)
            vel[vel == BAD_VEL] = np.nan
            vel /= 1000.0
            cor = amp = None
            if ids['cor'] in blocks:
                cor = np.frombuffer(b, 'u1', nc * 4, blocks[ids['cor']] + 2).reshape(nc, 4).astype(float)
            if ids['amp'] in blocks:
                amp = np.frombuffer(b, 'u1', nc * 4, blocks[ids['amp']] + 2).reshape(nc, 4).astype(float)
            vl = blocks.get(ids['var'])
            ens = struct.unpack_from('<H', b, vl + 2)[0] if vl else -1
            temp = struct.unpack_from('<h', b, vl + 26)[0] / 100.0 if vl else np.nan
            sspd = struct.unpack_from('<H', b, vl + 14)[0] if vl else np.nan
            bt = np.nan
            if 0x0600 in blocks:
                bo = blocks[0x0600]
                rng = np.array(struct.unpack_from('<4H', b, bo + 16), float) / 100.0
                rng[rng <= 0] = np.nan
                if np.isfinite(rng).sum() >= 3:
                    bt = np.nanmean(rng)
            t = dday.get(i, np.nan) if dday else np.nan
            recs.append((t, ens, fl, vel, cor, amp, temp, sspd, bt))
        i += nbytes + 2
    return recs


# ----------------------------------------------------------------------------- beam -> earth
def beam_to_inst(vb, beam_angle=30.0, three_beam=True):
    """vb (..., 4) beam velocities -> x, y, z, e (RDI convex janus)."""
    th = math.radians(beam_angle)
    a = 1.0 / (2.0 * math.sin(th)); bb = 1.0 / (4.0 * math.cos(th)); d = a / math.sqrt(2.0)
    b1, b2, b3, b4 = [vb[..., k].copy() for k in range(4)]
    if three_beam:
        # 3-beam solutions: set error velocity to zero and solve for the missing beam
        nbad = np.isnan(vb).sum(axis=-1)
        one = nbad == 1
        m = one & np.isnan(b1); b1[m] = (-b2 + b3 + b4)[m]
        m = one & np.isnan(b2); b2[m] = (-b1 + b3 + b4)[m]
        m = one & np.isnan(b3); b3[m] = (b1 + b2 - b4)[m]
        m = one & np.isnan(b4); b4[m] = (b1 + b2 - b3)[m]
    x = a * (b1 - b2)
    y = a * (b4 - b3)
    z = bb * (b1 + b2 + b3 + b4)
    e = d * (b1 + b2 - b3 - b4)
    return x, y, z, e


# ----------------------------------------------------------------------------- raw processing
DEFAULT_RAW_OPTS = dict(
    pingtype='bb',            # 'bb' or 'nb'
    inst='os75',
    enslen=300.0,             # s
    h_align=None,             # deg, None -> from proc.toml (h_align) or PD0 EA
    ducer_depth=None,         # m, None -> proc.toml
    heading_source='seapath', # 'seapath' | 'gyro' | 'gyro+hcorr'
    beam_angle=30.0,
    three_beam=True,
    min_corr=None,            # optional correlation threshold (counts), None = off
    max_err_vel=0.5,          # m/s single-ping error velocity limit (|e|), None = off
    max_w=0.5,                # m/s |w| single-ping limit, None = off
    pg_min=50.0,              # % good pings for an ensemble bin to be kept
    bottom_mask=True,         # mask bins below side-lobe contaminated depth (BT)
    xducer_dx=0.0, xducer_dy=0.0,   # transducer offset from GPS (m, stbd/fwd)
    sign=1.0,                 # sign convention of beam velocities (verified against CODAS)
    ens_edges=None,           # optional explicit ensemble edges (dday), e.g. to match CODAS
)


def process_raw(cf, opts=None, progress=None):
    o = dict(DEFAULT_RAW_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None or k in ('min_corr', 'max_err_vel', 'max_w', 'ens_edges')})
    inst = o['inst']
    pt = o['pingtype']
    log = [f'Raw single-ping processing: {inst}{pt}']
    files = cf.raw_files(inst)
    if not files:
        raise FileNotFoundError(f'No raw/{inst}/*.raw files')

    def prog(msg, frac):
        if progress:
            progress(msg, frac)

    # --- ADCP pings
    recs = []
    for k, f in enumerate(files):
        prog(f'PD0: {os.path.basename(f)}', 0.5 * k / len(files))
        recs.extend(read_pd0_file(f, pt))
    recs = [r for r in recs if np.isfinite(r[0])]
    if not recs:
        raise RuntimeError(f'No {pt} pings with timestamps found (check .raw.log files)')
    recs.sort(key=lambda r: r[0])
    t = np.array([r[0] for r in recs])
    maxnc = max(r[2]['ncell'] for r in recs)
    npg = len(recs)
    V = np.full((npg, maxnc, 4), np.nan)
    C = np.full((npg, maxnc, 4), np.nan)
    A = np.full((npg, maxnc, 4), np.nan)
    cfgkey = []
    for k, r in enumerate(recs):
        nc = r[2]['ncell']
        V[k, :nc] = r[3]
        if r[4] is not None:
            C[k, :nc] = r[4]
        if r[5] is not None:
            A[k, :nc] = r[5]
        cfgkey.append((nc, r[2]['cell'], r[2]['bin1'], r[2]['blank']))
    temp = np.array([r[6] for r in recs])
    bt_rng = np.array([r[8] for r in recs])
    ea = recs[0][2]['ea']
    log.append(f'  {npg} {pt} pings from {len(files)} files; configs: {sorted(set(cfgkey))}')

    # BT is only in the nb pings on an OS; borrow it if we are processing bb
    if pt == 'bb' and not np.isfinite(bt_rng).any():
        nbrec = []
        for f in files:
            nbrec.extend([(r[0], r[8]) for r in read_pd0_file(f, 'nb') if np.isfinite(r[0])])
        if nbrec:
            nbrec.sort()
            tb = np.array([x[0] for x in nbrec]); rb = np.array([x[1] for x in nbrec])
            ok = np.isfinite(rb)
            if ok.sum() > 1:
                bt_rng = np.interp(t, tb[ok], rb[ok], left=np.nan, right=np.nan)
                # don't interpolate across long gaps
                gap = np.interp(t, tb[ok], np.r_[0, np.diff(tb[ok])] * 86400, left=np.inf, right=np.inf)
                bt_rng[gap > 120] = np.nan
                log.append('  bottom track borrowed from interleaved nb pings')

    # --- heading
    prog('reading navigation', 0.55)
    sea = read_seapath(cf.nav_files('seapath', 'sea'))
    if o['heading_source'] in ('gyro', 'gyro+hcorr'):
        gy = read_gyro(cf.nav_files('gyro', 'hdg'))
        hdg = interp_angle(t, gy['t'], gy['hdg'])
        if o['heading_source'] == 'gyro+hcorr' and len(sea['th']):
            hs = interp_angle(gy['t'], sea['th'], sea['hdg'])
            dh = angdiff(hs, gy['hdg'])
            ok = np.isfinite(dh)
            # smooth correction over the ensemble length
            corr = _running_median_on_grid(gy['t'][ok], dh[ok], o['enslen'])
            hdg = (hdg + np.interp(t, gy['t'][ok], corr)) % 360
            log.append(f'  heading: gyro + seapath correction (median {np.nanmedian(dh):.2f} deg)')
        else:
            log.append('  heading: gyro')
    else:
        hdg = interp_angle(t, sea['th'], sea['hdg'])
        log.append('  heading: seapath INHDT')

    h_align = o['h_align']
    if h_align is None:
        h_align = cf.cfg_get('h_align', inst, None)
        if h_align is None:
            h_align = ea
    ducer = o['ducer_depth'] if o['ducer_depth'] is not None else float(cf.cfg_get('ducer_depth', inst, 5.0) or 5.0)
    log.append(f'  h_align = {h_align} deg, transducer depth = {ducer} m, beam angle = {o["beam_angle"]}')

    # --- beam -> instrument -> earth (ship-relative water velocity)
    Vb = V * o['sign']
    if o['min_corr'] is not None:
        Vb[C < o['min_corr']] = np.nan
    x, y, z, e = beam_to_inst(Vb, o['beam_angle'], o['three_beam'])
    H = np.radians((hdg + h_align) % 360)[:, None]
    ue = x * np.cos(H) + y * np.sin(H)
    ve = -x * np.sin(H) + y * np.cos(H)
    bad = ~np.isfinite(ue) | ~np.isfinite(hdg)[:, None]
    if o['max_err_vel'] is not None:
        bad |= np.abs(e) > o['max_err_vel']
    if o['max_w'] is not None:
        bad |= np.abs(z) > o['max_w']
    ue[bad] = np.nan; ve[bad] = np.nan; z = np.where(bad, np.nan, z); e = np.where(bad, np.nan, e)

    # --- positions & ship velocity at the transducer
    prog('ship velocity', 0.65)
    lat = np.interp(t, sea['t'], sea['lat'], left=np.nan, right=np.nan)
    lon = np.interp(t, sea['t'], sea['lon'], left=np.nan, right=np.nan)
    if o['xducer_dx'] or o['xducer_dy']:
        hs_ = np.radians(hdg)
        de = o['xducer_dx'] * np.cos(hs_) + o['xducer_dy'] * np.sin(hs_)
        dn = -o['xducer_dx'] * np.sin(hs_) + o['xducer_dy'] * np.cos(hs_)
        lat = lat + dn / 111194.9
        lon = lon + de / (111194.9 * np.cos(np.radians(lat)))
        log.append(f'  positions translated to transducer: dx={o["xducer_dx"]} dy={o["xducer_dy"]} m')

    # --- ensemble averaging
    prog('ensemble averaging', 0.75)
    enslen = o['enslen'] / 86400.0
    keys = np.array([hash(k) for k in cfgkey])
    if o.get('ens_edges') is not None:
        ens_id = np.searchsorted(np.asarray(o['ens_edges']), t)
    else:
        ens_id = np.floor((t - t[0]) / enslen).astype(int)
    brk = np.r_[True, (np.diff(ens_id) != 0) | (np.diff(keys) != 0)]
    grp = np.cumsum(brk) - 1
    ng = grp[-1] + 1
    out = {k: [] for k in ('time', 'lon', 'lat', 'u', 'v', 'w', 'e', 'amp', 'pg', 'uship', 'vship',
                           'heading', 'temp', 'np', 'depth', 'bt')}
    Ab = np.nanmean(A, axis=2)
    for g in range(ng):
        s = np.where(grp == g)[0]
        if len(s) < 2:
            continue
        ts = t[s]
        nc, cell, bin1, blank = cfgkey[s[0]]
        # ship velocity: displacement between first and last ping positions of the ensemble
        la, lo = lat[s], lon[s]
        okp = np.isfinite(la) & np.isfinite(lo)
        if okp.sum() >= 2:
            i0, i1 = s[okp][0], s[okp][-1]
            dtsec = (t[i1] - t[i0]) * 86400
            if dtsec > 0:
                us = (lon[i1] - lon[i0]) * 111194.9 * np.cos(np.radians(0.5 * (lat[i0] + lat[i1]))) / dtsec
                vs = (lat[i1] - lat[i0]) * 111194.9 / dtsec
            else:
                us = vs = np.nan
        else:
            us = vs = np.nan
        with np.errstate(all='ignore'):
            ngood = np.isfinite(ue[s]).sum(axis=0).astype(float)
            pg = 100.0 * ngood / len(s)
            ur = np.nanmean(ue[s], axis=0); vr = np.nanmean(ve[s], axis=0)
            out['w'].append(np.nanmean(z[s], axis=0)); out['e'].append(np.nanmean(e[s], axis=0))
            out['amp'].append(np.nanmean(Ab[s], axis=0))
        out['u'].append(ur + us); out['v'].append(vr + vs)
        out['pg'].append(pg); out['uship'].append(us); out['vship'].append(vs)
        out['time'].append(np.mean(ts))
        out['lon'].append(np.nanmean(lo) if okp.any() else np.nan)
        out['lat'].append(np.nanmean(la) if okp.any() else np.nan)
        hr = np.radians(hdg[s][np.isfinite(hdg[s])])
        out['heading'].append((np.degrees(np.arctan2(np.sin(hr).mean(), np.cos(hr).mean())) % 360) if len(hr) else np.nan)
        out['temp'].append(np.nanmean(temp[s])); out['np'].append(len(s))
        hh = hdg[s][np.isfinite(hdg[s])]
        out.setdefault('hrange', []).append(np.max(np.abs(angdiff(hh, hh[0]))) * 2 if len(hh) else np.nan)
        dep = ducer + bin1 + np.arange(maxnc) * cell
        dep[nc:] = np.nan
        out['depth'].append(dep)
        with np.errstate(all='ignore'):
            out['bt'].append(ducer + np.nanmedian(bt_rng[s]) * math.cos(math.radians(o['beam_angle']))
                             if np.isfinite(bt_rng[s]).any() else np.nan)
    arr = {k: np.array(v) for k, v in out.items()}
    hrange = arr.pop('hrange', None)
    sec = Section(name=f'{inst}{pt}', source='raw single-ping', yearbase=cf.config.get('yearbase', YEARBASE_DEFAULT),
                  time=arr['time'], lon=arr['lon'], lat=arr['lat'], depth=arr['depth'],
                  u=arr['u'], v=arr['v'], w=arr['w'], e=arr['e'], amp=arr['amp'], pg=arr['pg'],
                  uship=arr['uship'], vship=arr['vship'], heading=arr['heading'], tr_temp=arr['temp'],
                  num_pings=arr['np'], bt_depth=arr['bt'],
                  meta=dict(opts=o, h_align=h_align, ducer=ducer, heading_range=hrange), log=log)
    # PG and bottom editing are applied by apply_editing; set a pflag for bottom here
    if o['bottom_mask']:
        # side-lobe: reject bins deeper than ducer + (bottom - ducer)*cos(beam angle) - 1 bin
        cell = np.nanmedian(np.diff(sec.depth, axis=1))
        lim = ducer + (sec.bt_depth - ducer) * math.cos(math.radians(o['beam_angle'])) - 0.5 * cell
        sec.pflag = np.where(sec.depth > lim[:, None], 2, 0).astype(float)
        sec.log.append('  bottom: side-lobe mask flagged (pflag=2) from bottom track')
    sec.log.append(f'  {len(sec.time)} ensembles of {o["enslen"]:.0f} s')
    prog('done', 1.0)
    return sec


def _running_median_on_grid(t, x, win_s):
    """Median of x within +/- win/2 on the same grid (coarse, fast)."""
    if len(t) < 3:
        return x
    w = win_s / 86400.0
    edges = np.arange(t[0], t[-1] + w, w)
    idx = np.digitize(t, edges)
    med = np.full(len(edges) + 1, np.nan)
    for k in np.unique(idx):
        med[k] = np.nanmedian(x[idx == k])
    return med[idx]


# ----------------------------------------------------------------------------- calibration & editing
def apply_calibration(sec, rotate_deg=0.0, scale=1.0, log=True):
    """Rotate & scale ship-relative velocity; ocean = scale*R(rel) + ship.
    Uses uship/vship stored in the section."""
    out = sec.copy()
    us = np.nan_to_num(out.uship)[:, None]; vs = np.nan_to_num(out.vship)[:, None]
    ur = out.u - us; vr = out.v - vs
    ur, vr = rotate_uv(ur, vr, rotate_deg)
    out.u = scale * ur + us; out.v = scale * vr + vs
    if log:
        out.log.append(f'Calibration applied: rotate {rotate_deg:+.3f} deg, scale {scale:.4f}')
    return out


def apply_editing(sec, pg_min=50.0, use_pflag=True, zmin=None, zmax=None, tmin=None, tmax=None,
                  max_speed=None, spike_thresh=None, max_abs_e=None, bottom_frac=None, max_turn=None):
    """Recompute sec.mask from scratch. Returns a list of text lines describing what was removed."""
    m = ~np.isfinite(sec.u) | ~np.isfinite(sec.v)
    notes = [f'  missing/NaN: {m.sum()}']
    def add(cond, label):
        nonlocal m
        c = cond & ~m
        m = m | cond
        notes.append(f'  {label}: {c.sum()}')
    if pg_min is not None:
        add(np.nan_to_num(sec.pg, nan=0) < pg_min, f'PG < {pg_min:g}')
    if use_pflag:
        add(np.nan_to_num(sec.pflag) > 0, 'flagged (pflag)')
    if zmin is not None:
        add(sec.depth < zmin, f'depth < {zmin:g}')
    if zmax is not None:
        add(sec.depth > zmax, f'depth > {zmax:g}')
    if tmin is not None:
        add(np.repeat((sec.time < tmin)[:, None], sec.u.shape[1], 1), 'before tmin')
    if tmax is not None:
        add(np.repeat((sec.time > tmax)[:, None], sec.u.shape[1], 1), 'after tmax')
    if max_speed is not None:
        with np.errstate(invalid='ignore'):
            add(np.hypot(sec.u, sec.v) > max_speed, f'speed > {max_speed:g}')
    if max_abs_e is not None and np.isfinite(sec.e).any():
        with np.errstate(invalid='ignore'):
            add(np.abs(sec.e) > max_abs_e, f'|err vel| > {max_abs_e:g}')
    if bottom_frac is not None and np.isfinite(sec.bt_depth).any():
        add(sec.depth > (bottom_frac * sec.bt_depth)[:, None], f'deeper than {bottom_frac:g} x bottom')
    if max_turn is not None:
        hr = sec.meta.get('heading_range')
        if hr is None or len(hr) != len(sec.time):
            # CODAS: use heading change to neighbouring ensembles as a proxy
            h = sec.heading
            d = np.abs(angdiff(h[1:], h[:-1]))
            hr = np.fmax(np.r_[d, np.nan], np.r_[np.nan, d])
        add(np.repeat((np.nan_to_num(hr) > max_turn)[:, None], sec.u.shape[1], 1), f'heading change > {max_turn:g} deg (turns)')
    if spike_thresh is not None:
        # vertical 2nd-difference despike (CODAS-like "glitch" test)
        for comp, name in ((sec.u, 'u'), (sec.v, 'v')):
            x = np.where(m, np.nan, comp)
            d2 = np.full_like(x, np.nan)
            d2[:, 1:-1] = np.abs(x[:, 1:-1] - 0.5 * (x[:, :-2] + x[:, 2:]))
            with np.errstate(invalid='ignore'):
                add(d2 > spike_thresh, f'{name} vertical spike > {spike_thresh:g}')
    sec.mask = m
    notes.insert(0, f'Editing: {(~m).sum()} good of {m.size} cells')
    sec.log.extend(notes)
    return notes


def subset_time(sec, keep):
    """Return a copy of the section with only the profiles where keep is True."""
    out = sec.copy()
    keep = np.asarray(keep, bool)
    n = len(sec.time)
    for k, v in vars(out).items():
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n:
            setattr(out, k, v[keep])
    hr = out.meta.get('heading_range')
    if hr is not None and len(hr) == n:
        out.meta['heading_range'] = np.asarray(hr)[keep]
    return out


# ----------------------------------------------------------------------------- export
def export_csv(sec, path):
    dts = dday_to_datetime(sec.time, sec.yearbase)
    dist = sec.distance_km()
    with open(path, 'w') as f:
        f.write('datetime_utc,dday,lon,lat,distance_km,depth_m,u_ms,v_ms,speed_ms,dir_deg_to,w_ms,amp,pg,uship,vship,heading\n')
        spd, dr = uv_to_spd_dir(sec.um, sec.vm)
        for i in range(len(sec.time)):
            for j in range(sec.u.shape[1]):
                if sec.mask[i, j] or not np.isfinite(sec.depth[i, j]):
                    continue
                f.write(f'{dts[i]:%Y-%m-%dT%H:%M:%S},{sec.time[i]:.6f},{sec.lon[i]:.6f},{sec.lat[i]:.6f},{dist[i]:.3f},'
                        f'{sec.depth[i, j]:.2f},{sec.u[i, j]:.4f},{sec.v[i, j]:.4f},{spd[i, j]:.4f},{dr[i, j]:.1f},'
                        f'{sec.w[i, j]:.4f},{sec.amp[i, j]:.1f},{sec.pg[i, j]:.0f},{sec.uship[i]:.3f},'
                        f'{sec.vship[i]:.3f},{sec.heading[i]:.2f}\n')


def export_mat(sec, path):
    from scipy.io import savemat
    savemat(path, dict(time=sec.time, yearbase=sec.yearbase, lon=sec.lon, lat=sec.lat, depth=sec.depth,
                       u=sec.um, v=sec.vm, w=np.where(sec.mask, np.nan, sec.w), amp=sec.amp, pg=sec.pg,
                       uship=sec.uship, vship=sec.vship, heading=sec.heading, tr_temp=sec.tr_temp,
                       bt_depth=sec.bt_depth, mask=sec.mask.astype(np.uint8),
                       processing_log='\n'.join(sec.log)), do_compression=True)


def export_nc(sec, path):
    from scipy.io import netcdf_file
    nt, nz = sec.u.shape
    f = netcdf_file(path, 'w', version=2)
    f.title = f'{sec.name} processed velocities'
    f.source = sec.source
    f.history = f'Created {dt.datetime.utcnow():%Y-%m-%d %H:%M} UTC by Sruth VMADCP processor'
    f.processing_log = '\n'.join(sec.log)
    f.yearbase = int(sec.yearbase)
    f.createDimension('time', nt); f.createDimension('depth_cell', nz)
    FV = 1e38

    def var(name, data, dims, units, long):
        vv = f.createVariable(name, 'f8' if name == 'time' else 'f4', dims)
        d = np.array(data, float)
        d[~np.isfinite(d)] = FV
        vv[:] = d
        vv.units = units; vv.long_name = long; vv._FillValue = FV if name != 'time' else FV
    var('time', sec.time, ('time',), f'days since {sec.yearbase}-01-01 00:00:00', 'Decimal day')
    var('lon', sec.lon, ('time',), 'degrees_east', 'Longitude')
    var('lat', sec.lat, ('time',), 'degrees_north', 'Latitude')
    var('depth', sec.depth, ('time', 'depth_cell'), 'meter', 'Depth')
    var('u', sec.um, ('time', 'depth_cell'), 'meter second-1', 'Zonal velocity component')
    var('v', sec.vm, ('time', 'depth_cell'), 'meter second-1', 'Meridional velocity component')
    var('amp', sec.amp, ('time', 'depth_cell'), 'counts', 'Received signal strength')
    var('pg', sec.pg, ('time', 'depth_cell'), 'percent', 'Percent good pings')
    var('uship', sec.uship, ('time',), 'meter second-1', 'Ship zonal velocity')
    var('vship', sec.vship, ('time',), 'meter second-1', 'Ship meridional velocity')
    var('heading', sec.heading, ('time',), 'degrees', 'Ship heading')
    var('bt_depth', sec.bt_depth, ('time',), 'meter', 'Bottom depth from bottom track')
    f.close()
