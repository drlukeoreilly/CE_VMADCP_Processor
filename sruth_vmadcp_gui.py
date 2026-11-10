"""
Sruth VMADCP Processor - GUI for UHDAS / Teledyne Ocean Surveyor cruise folders.

    python sruth_vmadcp_gui.py [cruise_folder]

Point it at a UHDAS cruise folder (e.g. ...\\DMAND\\CE26011). It finds:
  proc/<sonar>/contour/<sonar>.nc     CODAS-processed ensembles (os75bb, os75nb)
  raw/os75/*.raw + .raw.log           single-ping PD0 data (bb + nb interleaved)
  raw/seapath/*.sea, raw/gyro/*.hdg   position / heading / attitude
  reports/ADCP_cals.txt               WT / BT / DXDY calibration estimates
and lets you calibrate, edit, plot and export.
Requires: numpy, scipy, matplotlib (tkinter ships with python.org Python).
"""
import os, sys, threading, traceback, datetime as dt
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import warnings
import numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.dates as mdates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sruth_vmadcp_core as core

APP = 'Sruth VMADCP Processor'
VARS = {'u (east)': 'u', 'v (north)': 'v', 'speed': 'spd', 'direction (to)': 'dir',
        'w (vertical)': 'w', 'amplitude': 'amp', 'percent good': 'pg', 'error velocity': 'e'}
XAXES = ['time', 'longitude', 'latitude', 'distance (km)', 'profile #']


class Scrollable(ttk.Frame):
    def __init__(self, parent, width=300):
        super().__init__(parent)
        c = tk.Canvas(self, width=width, highlightthickness=0)
        sb = ttk.Scrollbar(self, orient='vertical', command=c.yview)
        self.inner = ttk.Frame(c)
        self.inner.bind('<Configure>', lambda e: c.configure(scrollregion=c.bbox('all')))
        c.create_window((0, 0), window=self.inner, anchor='nw')
        c.configure(yscrollcommand=sb.set)
        c.pack(side='left', fill='both', expand=True); sb.pack(side='right', fill='y')
        c.bind_all('<MouseWheel>', lambda e: c.yview_scroll(int(-e.delta / 120), 'units') if self._over(e) else None)
        self._c = c

    def _over(self, e):
        w = e.widget
        while w is not None:
            if w == self:
                return True
            w = getattr(w, 'master', None)
        return False


class PlotTab(ttk.Frame):
    def __init__(self, parent, ctrl_builder=None):
        super().__init__(parent)
        self.ctrl = ttk.Frame(self); self.ctrl.pack(side='top', fill='x', padx=4, pady=2)
        self.fig = Figure(figsize=(9, 6), dpi=100, layout='constrained')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        tb = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False); tb.update()
        tb.pack(side='bottom', fill='x')
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)


class App(tk.Tk):
    def __init__(self, folder=None):
        super().__init__()
        self.title(APP)
        self.geometry('1450x900')
        self.cf = None
        self.base = None      # section as loaded (before cal/edit)
        self.sec = None       # current processed section
        self.compare = None   # optional second section
        self.nav = None
        self._build()
        if folder:
            self.folder_var.set(folder)
            self.after(200, self.scan_folder)

    # ------------------------------------------------------------------ layout
    def _build(self):
        top = ttk.Frame(self); top.pack(fill='x', padx=6, pady=4)
        ttk.Label(top, text='Cruise folder:').pack(side='left')
        self.folder_var = tk.StringVar()
        ttk.Entry(top, textvariable=self.folder_var, width=80).pack(side='left', fill='x', expand=True, padx=4)
        ttk.Button(top, text='Browse…', command=self.browse).pack(side='left')
        ttk.Button(top, text='Scan folder', command=self.scan_folder).pack(side='left', padx=4)

        pw = ttk.PanedWindow(self, orient='horizontal'); pw.pack(fill='both', expand=True)
        left = Scrollable(pw, width=330); pw.add(left, weight=0)
        right = ttk.Notebook(pw); pw.add(right, weight=1)
        self.nb = right
        L = left.inner

        # --- data source
        f = ttk.LabelFrame(L, text='1  Data source'); f.pack(fill='x', padx=4, pady=4)
        self.src_var = tk.StringVar()
        self.src_cb = ttk.Combobox(f, textvariable=self.src_var, state='readonly', width=38)
        self.src_cb.pack(fill='x', padx=4, pady=2)
        ttk.Button(f, text='Load / process', command=self.load_source).pack(fill='x', padx=4, pady=2)

        # --- raw options
        f = ttk.LabelFrame(L, text='2  Raw single-ping options'); f.pack(fill='x', padx=4, pady=4)
        self.r = {}
        def row(parent, label, key, default, width=10, kind='entry', values=None):
            fr = ttk.Frame(parent); fr.pack(fill='x', padx=4, pady=1)
            ttk.Label(fr, text=label, width=22).pack(side='left')
            if kind == 'entry':
                v = tk.StringVar(value=str(default)); ttk.Entry(fr, textvariable=v, width=width).pack(side='left')
            elif kind == 'check':
                v = tk.BooleanVar(value=default); ttk.Checkbutton(fr, variable=v).pack(side='left')
            else:
                v = tk.StringVar(value=default); ttk.Combobox(fr, textvariable=v, values=values, width=width, state='readonly').pack(side='left')
            return v
        self.r['enslen'] = row(f, 'Ensemble length (s)', 'enslen', 300)
        self.r['heading_source'] = row(f, 'Heading source', 'hs', 'gyro+hcorr', 12, 'combo', ['gyro+hcorr', 'seapath', 'gyro'])
        self.r['h_align'] = row(f, 'Transducer align (deg)', 'ha', '')
        self.r['ducer_depth'] = row(f, 'Transducer depth (m)', 'dd', '')
        self.r['xducer_dx'] = row(f, 'xducer dx (m, stbd)', 'dx', 0)
        self.r['xducer_dy'] = row(f, 'xducer dy (m, fwd)', 'dy', 0)
        ttk.Button(f, text='Use DXDY from ADCP_cals.txt', command=self.use_dxdy).pack(fill='x', padx=4, pady=1)
        self.r['max_err_vel'] = row(f, 'Ping |err vel| max (m/s)', 'ev', 0.5)
        self.r['max_w'] = row(f, 'Ping |w| max (m/s)', 'w', 0.5)
        self.r['three_beam'] = row(f, '3-beam solutions', 'tb', True, kind='check')
        self.r['bottom_mask'] = row(f, 'Side-lobe bottom flag', 'bm', True, kind='check')
        self.r['match_codas'] = row(f, 'Match CODAS ensembles', 'mc', False, kind='check')
        ttk.Label(f, text='Blank align/depth = from proc.toml', foreground='gray').pack(anchor='w', padx=4)

        # --- calibration
        f = ttk.LabelFrame(L, text='3  Calibration (ship-relative velocity)'); f.pack(fill='x', padx=4, pady=4)
        self.rot_var = tk.StringVar(value='0.0'); self.scale_var = tk.StringVar(value='1.0')
        fr = ttk.Frame(f); fr.pack(fill='x', padx=4)
        ttk.Label(fr, text='Rotate (deg)', width=14).pack(side='left'); ttk.Entry(fr, textvariable=self.rot_var, width=9).pack(side='left')
        fr = ttk.Frame(f); fr.pack(fill='x', padx=4)
        ttk.Label(fr, text='Scale factor', width=14).pack(side='left'); ttk.Entry(fr, textvariable=self.scale_var, width=9).pack(side='left')
        fr = ttk.Frame(f); fr.pack(fill='x', padx=4, pady=2)
        ttk.Button(fr, text='Use WT', command=lambda: self.use_cal('WT')).pack(side='left', expand=True, fill='x')
        ttk.Button(fr, text='Use BT', command=lambda: self.use_cal('BT')).pack(side='left', expand=True, fill='x')
        ttk.Button(fr, text='Reset', command=lambda: (self.rot_var.set('0.0'), self.scale_var.set('1.0'))).pack(side='left', expand=True, fill='x')
        self.cal_lbl = ttk.Label(f, text='', foreground='gray', wraplength=300, justify='left'); self.cal_lbl.pack(anchor='w', padx=4)

        # --- editing
        f = ttk.LabelFrame(L, text='4  Editing'); f.pack(fill='x', padx=4, pady=4)
        self.e = {}
        self.e['pg_min'] = row(f, 'PG minimum (%)', 'pg', 50)
        self.e['use_pflag'] = row(f, 'Apply CODAS/bottom flags', 'pf', True, kind='check')
        self.e['zmin'] = row(f, 'Depth min (m)', 'zmin', '')
        self.e['zmax'] = row(f, 'Depth max (m)', 'zmax', '')
        self.e['tmin'] = row(f, 'Start time (UTC)', 'tmin', '', width=17)
        self.e['tmax'] = row(f, 'End time (UTC)', 'tmax', '', width=17)
        self.trange_lbl = ttk.Label(f, text='Format: 2026-09-24 01:30  (or HH:MM, or dday)',
                                    foreground='gray', wraplength=300, justify='left')
        self.trange_lbl.pack(anchor='w', padx=4)
        fr = ttk.Frame(f); fr.pack(fill='x', padx=4)
        ttk.Button(fr, text='Fill full range', command=self.fill_trange).pack(side='left', expand=True, fill='x')
        ttk.Button(fr, text='Clear times', command=lambda: (self.e['tmin'].set(''), self.e['tmax'].set(''))).pack(side='left', expand=True, fill='x')
        self.e['max_speed'] = row(f, 'Max speed (m/s)', 'ms', 3.0)
        self.e['spike_thresh'] = row(f, 'Vertical spike (m/s)', 'sp', '')
        self.e['max_abs_e'] = row(f, 'Ens |err vel| max (m/s)', 'ae', '')
        self.e['bottom_frac'] = row(f, 'Keep < frac × BT depth', 'bf', '')
        self.e['max_turn'] = row(f, 'Reject turns > (deg/ens)', 'mt', '')
        ttk.Button(f, text='Apply calibration + editing', command=self.reprocess).pack(fill='x', padx=4, pady=4)

        # --- compare & export
        f = ttk.LabelFrame(L, text='5  Compare / export'); f.pack(fill='x', padx=4, pady=4)
        self.cmp_var = tk.StringVar()
        self.cmp_cb = ttk.Combobox(f, textvariable=self.cmp_var, state='readonly', width=38); self.cmp_cb.pack(fill='x', padx=4)
        ttk.Button(f, text='Load comparison dataset', command=self.load_compare).pack(fill='x', padx=4, pady=2)
        ttk.Separator(f).pack(fill='x', pady=3)
        self.x_csv = tk.BooleanVar(value=True); self.x_nc = tk.BooleanVar(value=True)
        self.x_mat = tk.BooleanVar(value=True); self.x_png = tk.BooleanVar(value=True)
        fr = ttk.Frame(f); fr.pack(fill='x', padx=4)
        for t, v in (('CSV', self.x_csv), ('NetCDF', self.x_nc), ('MAT', self.x_mat), ('PNGs', self.x_png)):
            ttk.Checkbutton(fr, text=t, variable=v).pack(side='left')
        ttk.Button(f, text='Export to <cruise>/sruth_processed', command=self.export).pack(fill='x', padx=4, pady=3)

        # progress
        self.prog = ttk.Progressbar(L, mode='determinate', maximum=1.0); self.prog.pack(fill='x', padx=6, pady=4)
        self.status = ttk.Label(L, text='Select a cruise folder', wraplength=310); self.status.pack(anchor='w', padx=6)

        # --- tabs
        self.tab_sum = ttk.Frame(right); right.add(self.tab_sum, text='Summary')
        self.sum_txt = tk.Text(self.tab_sum, wrap='none', font=('Consolas', 10)); self.sum_txt.pack(fill='both', expand=True)

        self.t_sec = PlotTab(right); right.add(self.t_sec, text='Section')
        c = self.t_sec.ctrl
        self.sec_vars = [tk.StringVar(value='u (east)'), tk.StringVar(value='v (north)')]
        ttk.Label(c, text='Panels:').pack(side='left')
        for v in self.sec_vars:
            ttk.Combobox(c, textvariable=v, values=['(none)'] + list(VARS), width=14, state='readonly').pack(side='left', padx=2)
        ttk.Label(c, text='  x:').pack(side='left')
        self.xaxis = tk.StringVar(value='time')
        ttk.Combobox(c, textvariable=self.xaxis, values=XAXES, width=12, state='readonly').pack(side='left')
        ttk.Label(c, text='  ±clim (m/s):').pack(side='left')
        self.clim = tk.StringVar(value='auto'); ttk.Entry(c, textvariable=self.clim, width=6).pack(side='left')
        self.show_masked = tk.BooleanVar(value=False)
        ttk.Checkbutton(c, text='show rejected', variable=self.show_masked).pack(side='left', padx=4)
        ttk.Button(c, text='Plot', command=self.plot_section).pack(side='left', padx=4)

        self.t_map = PlotTab(right); right.add(self.t_map, text='Map / vectors')
        c = self.t_map.ctrl
        ttk.Label(c, text='Layer (m):').pack(side='left')
        self.lz0 = tk.StringVar(value=''); self.lz1 = tk.StringVar(value='')
        ttk.Entry(c, textvariable=self.lz0, width=6).pack(side='left'); ttk.Label(c, text='to').pack(side='left')
        ttk.Entry(c, textvariable=self.lz1, width=6).pack(side='left')
        ttk.Label(c, text='  every n:').pack(side='left'); self.vsub = tk.StringVar(value='1')
        ttk.Entry(c, textvariable=self.vsub, width=4).pack(side='left')
        ttk.Label(c, text='  scale (deg per m/s):').pack(side='left'); self.vscale = tk.StringVar(value='auto')
        ttk.Entry(c, textvariable=self.vscale, width=6).pack(side='left')
        ttk.Button(c, text='Plot', command=self.plot_map).pack(side='left', padx=4)

        self.t_prof = PlotTab(right); right.add(self.t_prof, text='Profiles')
        c = self.t_prof.ctrl
        ttk.Label(c, text='Profile #').pack(side='left')
        self.prof_idx = tk.IntVar(value=0)
        self.prof_scale = ttk.Scale(c, from_=0, to=1, orient='horizontal', length=400,
                                    command=lambda v: (self.prof_idx.set(int(float(v))), self.plot_profiles()))
        self.prof_scale.pack(side='left', padx=4)
        self.prof_lbl = ttk.Label(c, text=''); self.prof_lbl.pack(side='left')

        self.t_ts = PlotTab(right); right.add(self.t_ts, text='Time series')
        c = self.t_ts.ctrl
        ttk.Label(c, text='Layer-mean current and ship data (uses Map layer)').pack(side='left')
        ttk.Button(c, text='Plot', command=self.plot_ts).pack(side='left', padx=4)

        self.t_nav = PlotTab(right); right.add(self.t_nav, text='Nav / heading QC')
        c = self.t_nav.ctrl
        ttk.Label(c, text='Reads raw/seapath & raw/gyro: heading difference, attitude, speed over ground').pack(side='left')
        ttk.Button(c, text='Load & plot', command=self.plot_nav).pack(side='left', padx=4)

        self.t_cmp = PlotTab(right); right.add(self.t_cmp, text='Compare')
        c = self.t_cmp.ctrl
        ttk.Label(c, text='Current vs comparison dataset (matched by nearest time & depth)').pack(side='left')
        ttk.Button(c, text='Plot', command=self.plot_compare).pack(side='left', padx=4)

        self.tab_log = ttk.Frame(right); right.add(self.tab_log, text='Log')
        self.log_txt = tk.Text(self.tab_log, wrap='none', font=('Consolas', 10)); self.log_txt.pack(fill='both', expand=True)

    # ------------------------------------------------------------------ helpers
    def set_status(self, msg, frac=None):
        def f():
            self.status.config(text=msg)
            if frac is not None:
                self.prog['value'] = frac
        self.after(0, f)

    def log(self, msg):
        def f():
            self.log_txt.insert('end', f'[{dt.datetime.now():%H:%M:%S}] {msg}\n'); self.log_txt.see('end')
        self.after(0, f)

    def run_bg(self, work, done):
        def target():
            try:
                res = work()
                self.after(0, lambda: done(res))
            except Exception as ex:
                tb = traceback.format_exc()
                self.log(tb)
                self.after(0, lambda: messagebox.showerror(APP, f'{ex}\n\nSee Log tab.'))
                self.set_status('Error', 0)
        threading.Thread(target=target, daemon=True).start()

    TIME_FORMATS = ('%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M',
                    '%d/%m/%Y %H:%M:%S', '%d/%m/%Y %H:%M', '%Y-%m-%d', '%d/%m/%Y')

    def _yearbase(self):
        if self.base is not None:
            return int(self.base.yearbase)
        return int(self.cf.config.get('yearbase', core.YEARBASE_DEFAULT)) if self.cf else core.YEARBASE_DEFAULT

    def tval(self, var, which='start'):
        """Parse a UTC time box -> decimal day. Accepts 'YYYY-MM-DD HH:MM[:SS]', 'DD/MM/YYYY HH:MM',
        a date alone, 'HH:MM' (on the dataset's first day, or next day if that is before the start), or dday."""
        s = str(var.get()).strip()
        if not s:
            return None
        yb = self._yearbase()
        base = dt.datetime(yb, 1, 1)
        try:
            return float(s)                         # plain decimal day still works
        except ValueError:
            pass
        for fmt in self.TIME_FORMATS:
            try:
                d = dt.datetime.strptime(s, fmt)
                if fmt in ('%Y-%m-%d', '%d/%m/%Y') and which == 'end':
                    d += dt.timedelta(days=1)       # a bare end date means end of that day
                return (d - base).total_seconds() / 86400.0
            except ValueError:
                continue
        for fmt in ('%H:%M:%S', '%H:%M'):
            try:
                hm = dt.datetime.strptime(s, fmt)
            except ValueError:
                continue
            if self.base is None:
                raise ValueError(f'"{s}": load a dataset first, or give a full date')
            t0 = self.base.time[0]
            day = np.floor(t0)
            dd = day + (hm.hour * 3600 + hm.minute * 60 + hm.second) / 86400.0
            if dd < t0 - 1e-6:
                dd += 1.0                           # time-of-day earlier than start -> next day
            return dd
        raise ValueError(f'time "{s}" not recognised - use e.g. 2026-09-24 01:30')

    def fill_trange(self):
        if self.base is None:
            messagebox.showinfo(APP, 'Load a dataset first'); return
        t = core.dday_to_datetime([self.base.time[0], self.base.time[-1]], self.base.yearbase)
        self.e['tmin'].set(f'{t[0]:%Y-%m-%d %H:%M:%S}'); self.e['tmax'].set(f'{t[1]:%Y-%m-%d %H:%M:%S}')

    def _update_trange_label(self):
        if self.base is None:
            return
        t = core.dday_to_datetime([self.base.time[0], self.base.time[-1]], self.base.yearbase)
        self.trange_lbl.config(text=f'Data: {t[0]:%Y-%m-%d %H:%M} to {t[1]:%Y-%m-%d %H:%M} UTC\n'
                                    f'Format: 2026-09-24 01:30  (or HH:MM, or dday)')

    @staticmethod
    def fval(var, default=None):
        s = str(var.get()).strip()
        if s == '' or s.lower() in ('none', 'auto'):
            return default
        return float(s)

    def browse(self):
        d = filedialog.askdirectory(title='Select UHDAS cruise folder')
        if d:
            self.folder_var.set(d); self.scan_folder()

    # ------------------------------------------------------------------ actions
    def scan_folder(self):
        try:
            self.cf = core.CruiseFolder(self.folder_var.get().strip())
        except Exception as ex:
            messagebox.showerror(APP, f'Cannot open folder:\n{ex}'); return
        opts = [f'{s}  (CODAS netCDF)' for s in self.cf.proc_sonars if self.cf.nc_file(s)]
        for inst in self.cf.raw_adcp_dirs:
            opts += [f'{inst}bb  (raw re-process)', f'{inst}nb  (raw re-process)']
        self.src_cb['values'] = opts; self.cmp_cb['values'] = opts
        if opts:
            self.src_var.set(opts[0])
            if len(opts) > 1:
                self.cmp_var.set(opts[1])
        self.sum_txt.delete('1.0', 'end'); self.sum_txt.insert('end', self.cf.summary() + '\n')
        for extra in ('uhdas_overview.txt', 'ADCP_settings.txt', 'ADCP_cals.txt'):
            p = os.path.join(self.cf.root, 'reports', extra)
            if os.path.exists(p):
                self.sum_txt.insert('end', f'\n===== reports/{extra} =====\n' + core._read_text(p))
        self.log(f'Scanned {self.cf.root}')
        self.set_status(f'Found {len(opts)} datasets. Choose one and press Load.', 0)
        self.update_cal_label()

    def _sonar_from(self, label):
        return label.split()[0]

    def _load(self, label):
        sonar = self._sonar_from(label)
        if '(CODAS' in label:
            self.set_status(f'Loading {sonar} CODAS…', 0.3)
            return core.load_codas(self.cf, sonar)
        inst, pt = sonar[:-2], sonar[-2:]
        o = dict(inst=inst, pingtype=pt,
                 enslen=self.fval(self.r['enslen'], 300.0),
                 heading_source=self.r['heading_source'].get(),
                 h_align=self.fval(self.r['h_align']), ducer_depth=self.fval(self.r['ducer_depth']),
                 xducer_dx=self.fval(self.r['xducer_dx'], 0.0), xducer_dy=self.fval(self.r['xducer_dy'], 0.0),
                 max_err_vel=self.fval(self.r['max_err_vel']), max_w=self.fval(self.r['max_w']),
                 three_beam=bool(self.r['three_beam'].get()), bottom_mask=bool(self.r['bottom_mask'].get()))
        if self.r['match_codas'].get() and self.cf.nc_file(sonar):
            tc = core.load_codas(self.cf, sonar).time          # CODAS time = ensemble end
            o['ens_edges'] = np.r_[tc[0] - o['enslen'] / 86400.0, tc]
        return core.process_raw(self.cf, o, progress=lambda m, fr: self.set_status(m, fr))

    def load_source(self):
        if not self.cf:
            self.scan_folder()
            if not self.cf:
                return
        label = self.src_var.get()
        if not label:
            return
        def done(sec):
            self.base = sec
            for line in sec.log:
                self.log(line)
            self.update_cal_label()
            self.reprocess()
            self.set_status(f'Loaded {sec.name} ({sec.source}): {len(sec.time)} profiles', 1.0)
        self.run_bg(lambda: self._load(label), done)

    def load_compare(self):
        label = self.cmp_var.get()
        if not (self.cf and label):
            return
        def done(sec):
            core.apply_editing(sec, pg_min=self.fval(self.e['pg_min'], 50.0), use_pflag=True)
            self.compare = sec
            self.log(f'Comparison loaded: {sec.name} ({sec.source})')
            self.nb.select(self.t_cmp); self.plot_compare()
            self.set_status('Comparison loaded', 1.0)
        self.run_bg(lambda: self._load(label), done)

    def current_sonar(self):
        return self.base.name if self.base else self._sonar_from(self.src_var.get() or 'os75bb')

    def update_cal_label(self):
        if not self.cf:
            return
        c = self.cf.cals.get(self.current_sonar(), {})
        if c:
            txt = '  '.join(f'{k}={v:g}' for k, v in c.items())
            self.cal_lbl.config(text=f'{self.current_sonar()} cals: {txt}\n(phase → Rotate, amp → Scale)')
        else:
            self.cal_lbl.config(text='No calibration estimates found for this sonar')

    def use_cal(self, kind):
        c = self.cf.cals.get(self.current_sonar(), {}) if self.cf else {}
        if f'{kind}_pha' not in c:
            messagebox.showinfo(APP, f'No {kind} calibration for {self.current_sonar()} in ADCP_cals.txt'); return
        self.rot_var.set(f'{c[kind + "_pha"]:.4f}'); self.scale_var.set(f'{c[kind + "_amp"]:.4f}')

    def use_dxdy(self):
        c = self.cf.cals.get(self._sonar_from(self.src_var.get() or 'os75bb'), {}) if self.cf else {}
        if 'dx' in c:
            self.r['xducer_dx'].set(f'{c["dx"]:.2f}'); self.r['xducer_dy'].set(f'{c["dy"]:.2f}')
        else:
            messagebox.showinfo(APP, 'No DXDY estimate found')

    def reprocess(self):
        if self.base is None:
            messagebox.showinfo(APP, 'Load a dataset first'); return
        try:
            tmin, tmax = self.tval(self.e['tmin'], 'start'), self.tval(self.e['tmax'], 'end')
            sec = core.apply_calibration(self.base, self.fval(self.rot_var, 0.0), self.fval(self.scale_var, 1.0))
            notes = core.apply_editing(sec, pg_min=self.fval(self.e['pg_min']), use_pflag=bool(self.e['use_pflag'].get()),
                                       zmin=self.fval(self.e['zmin']), zmax=self.fval(self.e['zmax']),
                                       tmin=tmin, tmax=tmax,
                                       max_speed=self.fval(self.e['max_speed']), spike_thresh=self.fval(self.e['spike_thresh']),
                                       max_abs_e=self.fval(self.e['max_abs_e']), bottom_frac=self.fval(self.e['bottom_frac']),
                                       max_turn=self.fval(self.e['max_turn']))
        except ValueError as ex:
            messagebox.showerror(APP, f'Check inputs: {ex}'); return
        self._update_trange_label()
        cal_line = sec.log[-len(notes) - 1]
        if tmin is not None or tmax is not None:
            fmt = lambda d: f'{core.dday_to_datetime(d, sec.yearbase)[0]:%Y-%m-%d %H:%M:%S}' if d is not None else 'open'
            keep = np.ones(len(sec.time), bool)
            if tmin is not None:
                keep &= sec.time >= tmin
            if tmax is not None:
                keep &= sec.time <= tmax
            if keep.sum() < 1:
                messagebox.showwarning(APP, 'No profiles inside that time window - check start/end times.'); return
            sec = core.subset_time(sec, keep)
            msg = f'Time clip: {fmt(tmin)} to {fmt(tmax)} UTC -> {keep.sum()} profiles kept'
            sec.log.append('  ' + msg); notes.append('  ' + msg)
        self.sec = sec
        self.prof_idx.set(0)
        self.log(cal_line)
        for n in notes:
            self.log(n)
        self.prof_scale.configure(to=max(len(sec.time) - 1, 1))
        self.sum_txt.insert('end', '\n===== current dataset =====\n' + sec.stats_text() + '\n' + '\n'.join(notes) + '\n')
        self.sum_txt.see('end')
        self.plot_section(); self.plot_map(); self.plot_profiles(); self.plot_ts()

    # ------------------------------------------------------------------ plotting
    def _x(self, sec):
        a = self.xaxis.get()
        if a == 'time':
            return mdates.date2num(core.dday_to_datetime(sec.time, sec.yearbase)), 'Time (UTC)', True
        if a == 'longitude':
            return sec.lon, 'Longitude (°E)', False
        if a == 'latitude':
            return sec.lat, 'Latitude (°N)', False
        if a == 'distance (km)':
            return sec.distance_km(), 'Distance along track (km)', False
        return np.arange(len(sec.time)), 'Profile #', False

    def _field(self, sec, key, masked=True):
        with np.errstate(all='ignore'):
            if key in ('u', 'v', 'w', 'e'):
                a = getattr(sec, key).copy()
            elif key == 'spd':
                a = np.hypot(sec.u, sec.v)
            elif key == 'dir':
                a = core.uv_to_spd_dir(sec.u, sec.v)[1]
            else:
                a = getattr(sec, key).astype(float).copy()
        if masked:
            a[sec.mask] = np.nan
        return a

    def plot_section(self):
        sec = self.sec
        fig = self.t_sec.fig; fig.clf()
        if sec is None:
            self.t_sec.canvas.draw(); return
        keys = [VARS[v.get()] for v in self.sec_vars if v.get() in VARS]
        if not keys:
            return
        x, xl, istime = self._x(sec)
        ok = np.isfinite(x)
        if self.xaxis.get() in ('longitude', 'latitude', 'distance (km)'):
            o = np.argsort(np.where(ok, x, np.inf));
        else:
            o = np.arange(len(x))
        axs = fig.subplots(len(keys), 1, sharex=True, squeeze=False)[:, 0]
        cl = self.fval(self.clim)
        for ax, k in zip(axs, keys):
            a = self._field(sec, k, masked=not self.show_masked.get())
            X = np.repeat(x[:, None], a.shape[1], 1)
            if k in ('u', 'v', 'w', 'e'):
                lim = cl if cl else max(np.nanpercentile(np.abs(a), 98), 0.05) if np.isfinite(a).any() else 1
                if k in ('w', 'e') and not cl:
                    lim = max(np.nanpercentile(np.abs(a), 98), 0.01) if np.isfinite(a).any() else 0.1
                kw = dict(cmap='RdBu_r', vmin=-lim, vmax=lim)
            elif k == 'spd':
                kw = dict(cmap='viridis', vmin=0, vmax=cl if cl else (np.nanpercentile(a, 98) if np.isfinite(a).any() else 1))
            elif k == 'dir':
                kw = dict(cmap='twilight', vmin=0, vmax=360)
            else:
                kw = dict(cmap='magma')
            sel = o[ok[o]]
            if len(sel) < 2:
                continue
            m = ax.pcolormesh(X[sel], self._filled_depth(sec)[sel], a[sel], shading='nearest', **kw)
            if np.isfinite(sec.bt_depth).any():
                ax.plot(x[sel], sec.bt_depth[sel], 'k-', lw=1.2, label='bottom (BT)')
            ax.set_ylabel('Depth (m)')
            ax.set_ylim(self._zbottom(sec, a), 0)
            units = {'u': 'm/s', 'v': 'm/s', 'w': 'm/s', 'e': 'm/s', 'spd': 'm/s', 'dir': '° (to)', 'amp': 'counts', 'pg': '%'}[k]
            fig.colorbar(m, ax=ax, label=units, pad=0.01)
            name = [n for n, kk in VARS.items() if kk == k][0]
            ax.set_title(f'{sec.name} — {name}  [{sec.source}]', fontsize=10, loc='left')
        if istime:
            axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
        axs[-1].set_xlabel(xl)
        self.t_sec.canvas.draw()

    @staticmethod
    def _zbottom(sec, a=None):
        """Plot depth limit: deepest finite value (+ margin), or bottom track if deeper."""
        a = sec.um if a is None else a
        z = np.where(np.isfinite(a), sec.depth, np.nan)
        zb = np.nanmax(z) if np.isfinite(z).any() else np.nanmax(sec.depth)
        if np.isfinite(sec.bt_depth).any():
            zb = max(zb, np.nanpercentile(sec.bt_depth, 95))
        return zb * 1.1 + 5

    @staticmethod
    def _filled_depth(sec):
        """pcolormesh needs finite coordinates: extend NaN depths beyond the last bin of each profile."""
        D = sec.depth.copy()
        for i in range(D.shape[0]):
            g = np.isfinite(D[i])
            if g.all():
                continue
            if g.sum() >= 2:
                idx = np.where(g)[0]
                dz = np.nanmedian(np.diff(D[i, idx]))
                k = np.arange(D.shape[1])
                D[i] = D[i, idx[0]] + (k - idx[0]) * dz
            else:
                D[i] = np.nanmedian(sec.depth, axis=0)
        return D

    def _layer(self, sec):
        z0 = self.fval(self.lz0); z1 = self.fval(self.lz1)
        if z0 is None:
            z0 = np.nanmin(sec.depth)
        if z1 is None:
            good = np.where(~sec.mask, sec.depth, np.nan)
            z1 = np.nanmax(good) if np.isfinite(good).any() else np.nanmax(sec.depth)
        return z0, z1

    def plot_map(self):
        sec = self.sec
        fig = self.t_map.fig; fig.clf()
        if sec is None:
            self.t_map.canvas.draw(); return
        ax = fig.add_subplot(111)
        z0, z1 = self._layer(sec)
        u, v = sec.layer_mean(z0, z1)
        n = max(int(self.fval(self.vsub, 1)), 1)
        ax.plot(sec.lon, sec.lat, '-', color='0.6', lw=1, zorder=1)
        ax.plot(sec.lon[0], sec.lat[0], 'go', ms=6, label='start'); ax.plot(sec.lon[-1], sec.lat[-1], 'rs', ms=6, label='end')
        s = slice(None, None, n)
        spd = np.hypot(u, v)
        sc = self.fval(self.vscale)
        span = max(np.nanmax(sec.lat) - np.nanmin(sec.lat), (np.nanmax(sec.lon) - np.nanmin(sec.lon)) * np.cos(np.radians(np.nanmean(sec.lat))), 0.01)
        if sc is None:
            ref = np.nanpercentile(spd, 90) if np.isfinite(spd).any() else 0.5
            sc = 0.15 * span / max(ref, 0.05)
        coslat = np.cos(np.radians(np.nanmean(sec.lat)))
        q = ax.quiver(sec.lon[s], sec.lat[s], u[s] * sc / coslat, v[s] * sc, spd[s], cmap='viridis',
                      angles='xy', scale_units='xy', scale=1, width=0.004, zorder=3)
        fig.colorbar(q, ax=ax, label='speed (m/s)')
        refv = 0.5
        ax.annotate(f'{refv} m/s', xy=(0.05, 0.04), xycoords='axes fraction')
        ax.set_aspect(1 / coslat, adjustable='datalim')
        ax.relim(); ax.autoscale_view()
        xl = ax.get_xlim(); yl = ax.get_ylim()
        ax.quiver([xl[0] + 0.05 * (xl[1] - xl[0])], [yl[0] + 0.1 * (yl[1] - yl[0])], [refv * sc / coslat], [0],
                  angles='xy', scale_units='xy', scale=1, width=0.004, color='k')
        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
        ax.set_title(f'{sec.name}: mean current {z0:.0f}–{z1:.0f} m', loc='left', fontsize=10)
        ax.legend(loc='upper right', fontsize=8); ax.grid(alpha=0.3)
        self.t_map.canvas.draw()

    def plot_profiles(self):
        sec = self.sec
        fig = self.t_prof.fig; fig.clf()
        if sec is None:
            self.t_prof.canvas.draw(); return
        i = min(max(self.prof_idx.get(), 0), len(sec.time) - 1)
        t = core.dday_to_datetime(sec.time[i], sec.yearbase)[0]
        self.prof_lbl.config(text=f'{i}  {t:%Y-%m-%d %H:%M:%S}  ({sec.lat[i]:.4f}, {sec.lon[i]:.4f})')
        z, mu, mv, su, sv = sec.mean_profile()
        a1, a2, a3 = fig.subplots(1, 3, sharey=True)
        for ax, m, s_, c, lab in ((a1, mu, su, 'C0', 'u'), (a1, mv, sv, 'C3', 'v')):
            ax.plot(m, z, c, lw=2, label=f'mean {lab}'); ax.fill_betweenx(z, m - s_, m + s_, color=c, alpha=0.15)
        a1.set_title('Mean ± std (all profiles)', fontsize=10); a1.set_xlabel('m/s'); a1.set_ylabel('Depth (m)'); a1.legend(fontsize=8)
        a2.plot(sec.um[i], sec.depth[i], 'C0.-', label='u'); a2.plot(sec.vm[i], sec.depth[i], 'C3.-', label='v')
        a2.plot(sec.u[i], sec.depth[i], 'C0', alpha=0.25); a2.plot(sec.v[i], sec.depth[i], 'C3', alpha=0.25)
        a2.set_title(f'Profile {i} (faint = rejected included)', fontsize=10); a2.set_xlabel('m/s'); a2.legend(fontsize=8)
        a3.plot(sec.pg[i], sec.depth[i], 'k.-', label='PG %')
        a3b = a3.twiny(); a3b.plot(sec.amp[i], sec.depth[i], 'C2.-', label='amp'); a3b.set_xlabel('amp', color='C2')
        a3.set_xlabel('Percent good'); a3.set_title('PG / amplitude', fontsize=10)
        if np.isfinite(sec.bt_depth[i]):
            for ax in (a1, a2, a3):
                ax.axhline(sec.bt_depth[i], color='k', ls='--', lw=1)
        for ax in (a1, a2, a3):
            ax.grid(alpha=0.3)
        a1.set_ylim(self._zbottom(sec, np.where(np.isfinite(sec.u), sec.u, np.nan)), 0)
        self.t_prof.canvas.draw()

    def plot_ts(self):
        sec = self.sec
        fig = self.t_ts.fig; fig.clf()
        if sec is None:
            self.t_ts.canvas.draw(); return
        z0, z1 = self._layer(sec)
        u, v = sec.layer_mean(z0, z1)
        tt = core.dday_to_datetime(sec.time, sec.yearbase)
        a1, a2, a3, a4 = fig.subplots(4, 1, sharex=True)
        a1.plot(tt, u, 'C0.-', label='u'); a1.plot(tt, v, 'C3.-', label='v'); a1.axhline(0, color='k', lw=0.5)
        a1.set_ylabel('m/s'); a1.set_title(f'{sec.name}: layer mean {z0:.0f}–{z1:.0f} m', loc='left', fontsize=10); a1.legend(fontsize=8)
        a2.plot(tt, np.hypot(u, v), 'k.-'); a2.set_ylabel('speed m/s')
        a2b = a2.twinx(); a2b.plot(tt, core.uv_to_spd_dir(u, v)[1], 'C1.', ms=4); a2b.set_ylabel('dir (to) °', color='C1'); a2b.set_ylim(0, 360)
        a3.plot(tt, np.hypot(sec.uship, sec.vship), 'k.-', label='ship speed'); a3.set_ylabel('ship m/s')
        a3b = a3.twinx(); a3b.plot(tt, sec.heading, 'C2.', ms=4); a3b.set_ylabel('heading °', color='C2'); a3b.set_ylim(0, 360)
        a4.plot(tt, sec.tr_temp, 'C4.-'); a4.set_ylabel('xducer T °C')
        if np.isfinite(sec.bt_depth).any():
            a4b = a4.twinx(); a4b.plot(tt, sec.bt_depth, 'k-', lw=1); a4b.set_ylabel('BT depth m'); a4b.invert_yaxis()
        a4.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
        for a in (a1, a2, a3, a4):
            a.grid(alpha=0.3)
        self.t_ts.canvas.draw()

    def plot_nav(self):
        if not self.cf:
            return
        def work():
            self.set_status('Reading seapath / gyro…', 0.2)
            sea = core.read_seapath(self.cf.nav_files('seapath', 'sea'), progress=lambda m, f: self.set_status(m, 0.2 + 0.6 * f))
            gy = core.read_gyro(self.cf.nav_files('gyro', 'hdg'))
            return sea, gy
        def done(res):
            sea, gy = res
            self.nav = res
            fig = self.t_nav.fig; fig.clf()
            yb = self.cf.config.get('yearbase', core.YEARBASE_DEFAULT)
            a1, a2, a3, a4 = fig.subplots(4, 1, sharex=True)
            if len(sea['th']):
                a1.plot(core.dday_to_datetime(sea['th'][::10], yb), sea['hdg'][::10], 'C0.', ms=1, label='seapath')
            if len(gy['t']):
                a1.plot(core.dday_to_datetime(gy['t'][::10], yb), gy['hdg'][::10], 'C1.', ms=1, label='gyro')
            a1.set_ylabel('heading °'); a1.legend(fontsize=8, markerscale=6)
            if len(gy['t']) and len(sea['th']):
                hs = core.interp_angle(gy['t'], sea['th'], sea['hdg'])
                dh = core.angdiff(hs, gy['hdg'])
                a2.plot(core.dday_to_datetime(gy['t'], yb), dh, 'k.', ms=1)
                med = np.nanmedian(dh)
                a2.axhline(med, color='r', lw=1)
                a2.set_ylim(med - 5, med + 5)
                a2.set_title(f'seapath − gyro heading (hcorr): median {med:.2f}°, std {np.nanstd(dh):.2f}°', fontsize=9, loc='left')
            a2.set_ylabel('Δheading °')
            if len(sea['tp']):
                tp = core.dday_to_datetime(sea['tp'][::5], yb)
                a3.plot(tp, sea['pitch'][::5], 'C0', lw=0.5, label='pitch'); a3.plot(tp, sea['roll'][::5], 'C3', lw=0.5, label='roll')
                a3.legend(fontsize=8); a3.set_ylabel('°')
            if len(sea['t']) > 2:
                t = sea['t']; la = sea['lat']; lo = sea['lon']
                dts = np.diff(t) * 86400
                sog = np.hypot(np.diff(lo) * 111194.9 * np.cos(np.radians(la[1:])), np.diff(la) * 111194.9) / np.where(dts > 0, dts, np.nan)
                k = max(len(sog) // 5000, 1)
                a4.plot(core.dday_to_datetime(t[1:][::k], yb), sog[::k], 'k', lw=0.5)
                a4.set_ylim(0, np.nanpercentile(sog, 99.5) * 1.2)
                gaps = np.where(dts > 10)[0]
                for g in gaps:
                    a4.axvspan(*core.dday_to_datetime([t[g], t[g + 1]], yb), color='r', alpha=0.2)
                a4.set_title(f'GPS speed over ground; {len(gaps)} gaps > 10 s shaded', fontsize=9, loc='left')
            a4.set_ylabel('SOG m/s')
            a4.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
            for a in (a1, a2, a3, a4):
                a.grid(alpha=0.3)
            self.t_nav.canvas.draw()
            self.set_status('Nav QC plotted', 1.0)
        self.nb.select(self.t_nav)
        self.run_bg(work, done)

    def plot_compare(self):
        A, B = self.sec, self.compare
        fig = self.t_cmp.fig; fig.clf()
        if A is None or B is None:
            ax = fig.add_subplot(111); ax.text(0.5, 0.5, 'Load a dataset and a comparison dataset', ha='center'); self.t_cmp.canvas.draw(); return
        # match B to A by nearest time (within 5 min) and interpolate in depth
        ta, tb = A.mid_time(), B.mid_time()
        ib = np.array([np.nanargmin(np.abs(tb - t)) for t in ta])
        okt = np.abs(tb[ib] - ta) < 120 / 86400
        Bu = np.full_like(A.u, np.nan); Bv = np.full_like(A.v, np.nan)
        for i in np.where(okt)[0]:
            j = ib[i]; zb = B.depth[j]; g = np.isfinite(zb) & np.isfinite(B.um[j])
            if g.sum() > 1:
                Bu[i] = np.interp(A.depth[i], zb[g], B.um[j][g], left=np.nan, right=np.nan)
                Bv[i] = np.interp(A.depth[i], zb[g], B.vm[j][g], left=np.nan, right=np.nan)
        du = A.um - Bu; dv = A.vm - Bv
        ok = np.isfinite(du) & np.isfinite(dv)
        a1, a2 = fig.subplots(1, 2)
        a1.plot(Bu[ok], A.um[ok], 'C0.', ms=3, label='u'); a1.plot(Bv[ok], A.vm[ok], 'C3.', ms=3, label='v')
        lim = np.nanmax(np.abs(np.r_[A.um[ok], Bu[ok], 0.1]))
        a1.plot([-lim, lim], [-lim, lim], 'k-', lw=0.7); a1.set_aspect('equal')
        a1.set_xlabel(f'{B.name} [{B.source}] m/s'); a1.set_ylabel(f'{A.name} [{A.source}] m/s'); a1.legend(); a1.grid(alpha=0.3)
        with np.errstate(all='ignore'):
            z = np.nanmedian(A.depth, 0)
            a2.plot(np.nanmean(du, 0), z, 'C0.-', label='mean Δu'); a2.plot(np.nanmean(dv, 0), z, 'C3.-', label='mean Δv')
            a2.plot(np.nanstd(du, 0), z, 'C0:', label='std Δu'); a2.plot(np.nanstd(dv, 0), z, 'C3:', label='std Δv')
        a2.axvline(0, color='k', lw=0.5); a2.invert_yaxis(); a2.set_xlabel('m/s'); a2.set_ylabel('Depth (m)'); a2.legend(fontsize=8); a2.grid(alpha=0.3)
        if ok.any():
            txt = (f'N={ok.sum()}  bias Δu={np.mean(du[ok]):+.3f}  Δv={np.mean(dv[ok]):+.3f}  '
                   f'rms Δu={np.sqrt(np.mean(du[ok] ** 2)):.3f}  Δv={np.sqrt(np.mean(dv[ok] ** 2)):.3f} m/s')
            fig.suptitle(txt, fontsize=10); self.log('Compare: ' + txt)
        self.t_cmp.canvas.draw()

    # ------------------------------------------------------------------ export
    def export(self):
        if self.sec is None:
            messagebox.showinfo(APP, 'Nothing to export'); return
        sec = self.sec
        tag = 'codas' if 'CODAS' in sec.source else 'raw'
        outdir = os.path.join(self.cf.root, 'sruth_processed', f'{sec.name}_{tag}_{dt.datetime.now():%Y%m%d_%H%M%S}')
        os.makedirs(outdir, exist_ok=True)
        base = os.path.join(outdir, f'{self.cf.config.get("cruiseid", self.cf.name)}_{sec.name}')
        written = []
        try:
            if self.x_csv.get():
                core.export_csv(sec, base + '.csv'); written.append(base + '.csv')
            if self.x_nc.get():
                core.export_nc(sec, base + '.nc'); written.append(base + '.nc')
            if self.x_mat.get():
                core.export_mat(sec, base + '.mat'); written.append(base + '.mat')
            if self.x_png.get():
                for name, tab in (('section', self.t_sec), ('map', self.t_map), ('profiles', self.t_prof),
                                  ('timeseries', self.t_ts), ('navqc', self.t_nav), ('compare', self.t_cmp)):
                    if tab.fig.axes:
                        p = f'{base}_{name}.png'; tab.fig.savefig(p, dpi=150); written.append(p)
            rep = base + '_processing_report.txt'
            with open(rep, 'w') as f:
                f.write(f'{APP}\nExported {dt.datetime.now():%Y-%m-%d %H:%M:%S}\nCruise folder: {self.cf.root}\n\n')
                f.write(sec.stats_text() + '\n\nSettings\n')
                f.write(f'  rotate = {self.rot_var.get()} deg, scale = {self.scale_var.get()}\n')
                for k, v in self.e.items():
                    f.write(f'  edit.{k} = {v.get()}\n')
                if 'raw' in sec.source:
                    for k, v in self.r.items():
                        f.write(f'  raw.{k} = {v.get()}\n')
                f.write('\nProcessing log\n' + '\n'.join(sec.log) + '\n')
            written.append(rep)
        except Exception as ex:
            self.log(traceback.format_exc()); messagebox.showerror(APP, f'Export failed: {ex}'); return
        self.log(f'Exported {len(written)} files to {outdir}')
        messagebox.showinfo(APP, f'Exported {len(written)} files to:\n{outdir}')


if __name__ == '__main__':
    App(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()
