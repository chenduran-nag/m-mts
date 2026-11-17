"""M-MTS: Adaptive Multi-Modal Trajectory Synthesis with Morphing Safety Envelopes.

Reference implementation used for the experiments in the paper, together with the
two baselines (FE-CBF and STCS).  Pure Python + NumPy.

Conventions
-----------
State of a vehicle  x = [px, py, psi, v]  (rear-axle position, heading, speed).
Inputs              u = [a, kappa]        (acceleration, path curvature).
Dynamics            px' = v cos psi, py' = v sin psi, psi' = v kappa, v' = a.
The footprint (4.0 m x 1.8 m) and the safety envelope are centred L/2 ahead of the
rear axle.
"""
import math, sys, time, heapq
import numpy as np

sys.setrecursionlimit(20000)

# ----------------------------------------------------------------------------
# Parameters (all values reported in the paper's parameter table)
# ----------------------------------------------------------------------------
P = dict(
    dt=0.1,            # control / simulation step [s]
    substeps=2,        # Euler sub-steps per control step
    plan_period=0.5,   # receding-horizon re-planning period [s]
    T=30.0,            # planning horizon [s] (covers the remaining trip)
    T_unc=6.0,         # prediction horizon for non-connected vehicles [s]
    dt_chk=0.25,       # conflict-check / profile sampling interval [s]
    L=2.5,             # wheelbase [m]
    car_len=4.0, car_wid=1.8,
    r=2.4,             # nominal envelope radius [m]
    r0=1.0,            # minimum lateral semi-axis [m]
    k_par=0.15,        # speed gain of longitudinal semi-axis [s/m]
    n_min=2.0, n_max=6.0,
    beta=10.0,         # sharpness of smooth min/max (log-sum-exp)
    lat_window=6.0,    # extent beyond the semi-axis within which neighbours shape d_lat [m]
    delta_ctx=0.45,    # clearance reserved when contracting b [m]
    gamma=1.0,         # CBF class-K gain
    a_min=-6.0, a_max=3.0, kappa_max=0.2,
    w_kappa=1000.0,    # QP weight on curvature deviation (braking preferred over swerving)
    sense=20.0,        # CBF neighbourhood radius [m]
    Delta=0.5, Dmax=10.0,  # DFCS delay step and maximum delay [s]
    node_cap=300,      # DFCS expansion budget per planning cycle
    sigma=0.5,         # mode-likelihood std [m]
    w_window=0.5,      # observation window for the mode-weight update [s]
    w_floor=0.02,      # lower bound on mode weights
    w_active=0.10,     # hypotheses with weight above this are checked in DFCS
    stcs_R=2.4,        # STCS sphere radius [m]
    stcs_delta=0.25,
    aging=0.2,         # priority gained per second spent yielding (one rank per 5 s)
    cost=dict(goal=1.0, time=0.5, effort=0.02, jerk=0.002, clear=10.0, c_ref=3.0),
)
MODES = [  # behavioural modes: (name, cruise speed [m/s], acceleration [m/s^2], margin delta [m])
    dict(name='cautious', v=4.0, acc=1.5, delta=0.35),
    dict(name='nominal', v=6.0, acc=2.0, delta=0.25),
    dict(name='aggressive', v=7.0, acc=3.0, delta=0.15),
]
NOMINAL = 1
HALF_L = P['car_len'] / 2
HALF_W = P['car_wid'] / 2
OFF = P['L'] / 2          # rear axle -> footprint centre


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------
M_VERT = 20
_ang = np.linspace(0, 2 * np.pi, M_VERT, endpoint=False)
_ca, _sa = np.cos(_ang), np.sin(_ang)


def superellipse(cx, cy, psi, a, b, n):
    e = 2.0 / n
    x = a * np.sign(_ca) * np.abs(_ca) ** e
    y = b * np.sign(_sa) * np.abs(_sa) ** e
    c, s = math.cos(psi), math.sin(psi)
    return np.stack([cx + c * x - s * y, cy + s * x + c * y], 1)


_rect = np.array([[1, -1], [1, 1], [-1, 1], [-1, -1]], float)


def rect(cx, cy, psi, hl=HALF_L, hw=HALF_W):
    x, y = _rect[:, 0] * hl, _rect[:, 1] * hw
    c, s = math.cos(psi), math.sin(psi)
    return np.stack([cx + c * x - s * y, cy + s * x + c * y], 1)


def _vert_edge(A, B):
    E = np.roll(B, -1, 0) - B
    d = A[:, None, :] - B[None, :, :]
    t = np.clip((d * E).sum(-1) / (E * E).sum(-1), 0, 1)
    pr = d - t[..., None] * E
    return math.sqrt((pr * pr).sum(-1).min())


def poly_sd(Pa, Qa):
    """Signed distance between two convex CCW polygons (negative = penetration depth)."""
    sep = -1e9
    for A, B in ((Pa, Qa), (Qa, Pa)):
        E = np.roll(A, -1, 0) - A
        N = np.stack([E[:, 1], -E[:, 0]], 1)
        N /= np.sqrt((N * N).sum(1))[:, None]
        s = (B @ N.T).min(0) - (A * N).sum(1)
        sep = max(sep, s.max())
    if sep <= 0:
        return sep
    return min(_vert_edge(Pa, Qa), _vert_edge(Qa, Pa))


def point_poly_sd(Pa, p):
    E = np.roll(Pa, -1, 0) - Pa
    d = p[None, :] - Pa
    t = np.clip((d * E).sum(1) / (E * E).sum(1), 0, 1)
    pr = d - t[:, None] * E
    dist = math.sqrt((pr * pr).sum(1).min())
    inside = np.all(E[:, 0] * d[:, 1] - E[:, 1] * d[:, 0] >= 0)
    return -dist if inside else dist


def poly_disc_sd(Pa, c, rad):
    return point_poly_sd(Pa, np.asarray(c, float)) - rad


def smin_arr(x, axis=-1, beta=P['beta']):
    m = x.min(axis=axis, keepdims=True)
    return (m - np.log(np.exp(-beta * (x - m)).sum(axis=axis, keepdims=True)) / beta).squeeze(axis)


def smin2(a, b, beta=P['beta']):
    return smin_arr(np.stack([a, b], -1), -1, beta)


def smax2(a, b, beta=P['beta']):
    return -smin2(-a, -b, beta)


# ----------------------------------------------------------------------------
# Morphing envelope (Eq. 7) -- vectorised over time samples
# ----------------------------------------------------------------------------
def _window_pen(xi, a):
    W = P['lat_window']
    return 3.0 * np.maximum(0, np.abs(xi) - (a + W)) ** 2


def dlat_static(cx, cy, psi, v, obs):
    """Lateral free space to static/known discs. cx..v: (T,), obs: (O,3) -> (T,)."""
    T = len(cx)
    big = np.full(T, P['r'] + 2.0)
    if obs is None or len(obs) == 0:
        return big
    a = P['r'] * (1 + P['k_par'] * np.maximum(v, 0))
    dx = obs[None, :, 0] - cx[:, None]
    dy = obs[None, :, 1] - cy[:, None]
    c, s = np.cos(psi)[:, None], np.sin(psi)[:, None]
    xi = c * dx + s * dy
    eta = -s * dx + c * dy
    cl = np.abs(eta) - obs[None, :, 2] - P['delta_ctx'] + _window_pen(xi, a[:, None])
    return smin_arr(np.concatenate([cl, big[:, None]], 1), 1)


def dlat_partner(cx, cy, psi, v, ox, oy, shared=True, orad=0.0):
    a = P['r'] * (1 + P['k_par'] * np.maximum(v, 0))
    dx, dy = ox - cx, oy - cy
    c, s = np.cos(psi), np.sin(psi)
    xi = c * dx + s * dy
    eta = -s * dx + c * dy
    if shared:
        cl = (np.abs(eta) - P['delta_ctx']) / 2
    else:
        cl = np.abs(eta) - orad - P['delta_ctx']
    return cl + _window_pen(xi, a)


def env_params(method, v, dlat):
    v = np.asarray(v, float)
    if method == 'fecbf':
        o = np.ones_like(v)
        return P['r'] * o, P['r'] * o, 2.0 * o
    a = P['r'] * (1 + P['k_par'] * np.maximum(v, 0))
    b = smax2(np.full_like(v, P['r0']), smin2(np.full_like(v, P['r']), dlat))
    n = P['n_min'] + (P['n_max'] - P['n_min']) * np.clip((P['r'] - b) / (P['r'] - P['r0']), 0, 1)
    return a, b, n


def centre(x):
    return x[0] + OFF * math.cos(x[2]), x[1] + OFF * math.sin(x[2])


# ----------------------------------------------------------------------------
# Paths (grid A* + shortcut smoothing) and speed profiles
# ----------------------------------------------------------------------------
LAT = dict(res=0.5, nh=16, kappa=0.15, step=1.0, clr=0.8, turn_cost=0.3)


def _foot_clear(x, y, h, obs):
    """Clearance between the footprint (rear axle at x, y, heading h) and the disc obstacles."""
    c, s_ = math.cos(h), math.sin(h)
    cx, cy = x + OFF * c, y + OFF * s_
    dx, dy = obs[:, 0] - cx, obs[:, 1] - cy
    xi, eta = np.abs(c * dx + s_ * dy), np.abs(-s_ * dx + c * dy)
    return (np.hypot(np.maximum(xi - HALF_L, 0), np.maximum(eta - HALF_W, 0)) - obs[:, 2]).min()


def _primitives():
    """Motion primitives of the state lattice: straight (1 m) and arcs of curvature +-kappa
    that change the heading by one heading bin (2*pi/nh)."""
    dpsi = 2 * math.pi / LAT['nh']
    L_arc = dpsi / LAT['kappa']
    prims = [(0.0, LAT['step'], 0)]
    for sgn in (1, -1):
        prims.append((sgn * LAT['kappa'], L_arc, sgn))
    out = []
    for k, L, dh in prims:
        ss = np.linspace(0, L, max(2, int(L / 0.25) + 1))
        if k == 0:
            pts = np.stack([ss, np.zeros_like(ss), np.zeros_like(ss)], 1)
        else:
            pts = np.stack([np.sin(k * ss) / k, (1 - np.cos(k * ss)) / k, k * ss], 1)
        out.append((pts, L + (LAT['turn_cost'] if k else 0.0), dh))
    return out


_PRIMS = _primitives()


def lattice_path(start, goal, obs, bounds):
    """State-lattice A* (x, y on a 0.5 m grid, 16 headings) with curvature-bounded primitives and
    footprint collision checking.  Returns a Path for the rear axle, or None."""
    if obs is None or len(obs) == 0:
        return straight(start, goal)
    x0, x1, y0, y1 = bounds
    res, nh = LAT['res'], LAT['nh']
    dpsi = 2 * math.pi / nh
    h0 = math.atan2(goal[1] - start[1], goal[0] - start[0])
    ih0 = int(round(h0 / dpsi)) % nh
    sx, sy = float(start[0]), float(start[1])
    key = lambda x, y, ih: (int(round(x / res)), int(round(y / res)), ih)
    hfun = lambda x, y: math.hypot(goal[0] - x, goal[1] - y)
    openq = [(hfun(sx, sy), 0.0, sx, sy, ih0)]
    came = {key(sx, sy, ih0): None}
    gbest = {key(sx, sy, ih0): 0.0}
    node = {key(sx, sy, ih0): (sx, sy, ih0, None)}
    found = None
    it = 0
    while openq and it < 200000:
        it += 1
        f, g, x, y, ih = heapq.heappop(openq)
        k0 = key(x, y, ih)
        if g > gbest[k0] + 1e-9:
            continue
        if hfun(x, y) < 1.5:
            found = k0
            break
        h = ih * dpsi
        c, s_ = math.cos(h), math.sin(h)
        for pts, cost, dh in _PRIMS:
            px = x + c * pts[:, 0] - s_ * pts[:, 1]
            py = y + s_ * pts[:, 0] + c * pts[:, 1]
            nx_, ny_ = px[-1], py[-1]
            if not (x0 <= nx_ <= x1 and y0 <= ny_ <= y1):
                continue
            ok = True
            for q in range(1, len(pts)):
                if _foot_clear(px[q], py[q], h + pts[q, 2], obs) < LAT['clr']:
                    ok = False
                    break
            if not ok:
                continue
            nih = (ih + dh) % nh
            k1 = key(nx_, ny_, nih)
            ng = g + cost
            if ng < gbest.get(k1, 1e18):
                gbest[k1] = ng
                came[k1] = k0
                node[k1] = (nx_, ny_, nih, (px, py))
                heapq.heappush(openq, (ng + hfun(nx_, ny_), ng, nx_, ny_, nih))
    if found is None:
        return None
    segs = []
    k = found
    while k is not None:
        seg = node[k][3]
        if seg is not None:
            segs.append(np.stack(seg, 1))
        k = came[k]
    segs.reverse()
    pts = np.concatenate([np.array([[sx, sy]])] + [sg[1:] for sg in segs] + [np.array([goal], float)], 0)
    return resample(pts)


astar_path = lattice_path


def straight(start, goal):
    return resample(np.array([start, goal], float))


def resample(pts, ds=0.5, smooth=False):
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)])
    n = max(2, int(s[-1] / ds) + 1)
    ss = np.linspace(0, s[-1], n)
    xy = np.stack([np.interp(ss, s, pts[:, 0]), np.interp(ss, s, pts[:, 1])], 1)
    if smooth and len(xy) > 7:
        for _ in range(3):
            xy[1:-1] = 0.25 * xy[:-2] + 0.5 * xy[1:-1] + 0.25 * xy[2:]
    return Path(xy)


class Path:
    def __init__(self, xy):
        self.xy = xy
        seg = np.hypot(*np.diff(xy, axis=0).T)
        self.s = np.concatenate([[0], np.cumsum(seg)])
        self.S = self.s[-1]
        d = np.diff(xy, axis=0)
        h = np.arctan2(d[:, 1], d[:, 0])
        self.h = np.unwrap(np.concatenate([h, h[-1:]]))

    def at(self, s):
        s = np.clip(s, 0, self.S)
        return (np.interp(s, self.s, self.xy[:, 0]), np.interp(s, self.s, self.xy[:, 1]),
                np.interp(s, self.s, self.h))

    def project(self, p, hint):
        lo, hi = max(0, hint - 30), min(len(self.xy) - 1, hint + 30)
        if np.min(np.hypot(*(self.xy[lo:hi + 1] - p).T)) > 3.0:   # far off: global search
            lo, hi = 0, len(self.xy) - 1
        best, bi, bt = 1e18, hint, 0.0
        for i in range(lo, max(lo + 1, hi)):
            a, b = self.xy[i], self.xy[min(i + 1, len(self.xy) - 1)]
            e = b - a
            ee = e @ e
            t = 0.0 if ee == 0 else min(1, max(0, (p - a) @ e / ee))
            q = a + t * e
            d = (p - q) @ (p - q)
            if d < best:
                best, bi, bt = d, i, t
        s = self.s[bi] + bt * (self.s[min(bi + 1, len(self.s) - 1)] - self.s[bi])
        return s, math.sqrt(best), bi


def speed_profile(s0, v0, S, mode, hold):
    """1-D profile from (s0, v0): hold (brake to rest and wait) for `hold` s, then cruise."""
    dt, n = P['dt_chk'], int(round(P['T'] / P['dt_chk'])) + 1
    s = np.empty(n)
    v = np.empty(n)
    s[0], v[0] = s0, v0
    for k in range(1, n):
        t = (k - 1) * dt
        vt = 0.0 if t < hold else mode['v']
        vt = min(vt, math.sqrt(max(0.0, 2 * 2.5 * (S - s[k - 1] - 0.3))))
        dv = min(max(vt - v[k - 1], -4.0 * dt), mode['acc'] * dt)
        v[k] = max(0.0, v[k - 1] + dv)
        s[k] = min(S, s[k - 1] + 0.5 * (v[k - 1] + v[k]) * dt)
    return s, v


# ----------------------------------------------------------------------------
# World objects
# ----------------------------------------------------------------------------
class Agent:
    def __init__(self, i, path, v0, prio):
        self.i, self.path, self.prio = i, path, prio
        x, y, h = path.at(0.0)
        self.x = np.array([float(x), float(y), float(h), v0])
        self.hint = 0
        self.s = 0.0
        self.done = False
        self.t_done = None
        self.plan = None      # dict(t0, s, v, mode)
        self.mode = NOMINAL
        self.track_err = []


class NCV:
    """Non-connected vehicle following a scripted trajectory fn(t) -> (x, y, psi, v)."""
    def __init__(self, fn, lane_psi):
        self.fn, self.lane_psi = fn, lane_psi
        self.w = np.ones(3) / 3

    def state(self, t):
        return np.array(self.fn(t), float)


class DynObs:
    def __init__(self, x0, y0, vx, vy, rad, t0=0.0):
        self.p0, self.vel, self.rad, self.t0 = np.array([x0, y0], float), np.array([vx, vy], float), rad, t0

    def pos(self, t):
        return self.p0 + self.vel * max(0.0, t - self.t0)

    def velocity(self, t):
        return self.vel if t >= self.t0 else np.zeros(2)


# ----------------------------------------------------------------------------
# Trajectory samples used by DFCS
# ----------------------------------------------------------------------------
class Traj:
    __slots__ = ('cx', 'cy', 'psi', 'v', 'dls', 'kind', 'rad', 's', 'vprof', 'delta', 'tag')

    def __init__(self, cx, cy, psi, v, dls, kind, rad=0.0, s=None, vprof=None, delta=0.25, tag=None):
        self.cx, self.cy, self.psi, self.v, self.dls = cx, cy, psi, v, dls
        self.kind, self.rad, self.s, self.vprof, self.delta, self.tag = kind, rad, s, vprof, delta, tag


def agent_traj(ag, mode_i, d_i, static_obs, method):
    mode = MODES[mode_i]
    s, v = speed_profile(ag.s, ag.x[3], ag.path.S, mode, d_i * P['Delta'])
    x, y, h = ag.path.at(s)
    cx, cy = x + OFF * np.cos(h), y + OFF * np.sin(h)
    gone = s >= ag.path.S - 1.0           # agents leave the scene on arrival
    cx, cy = np.where(gone, np.nan, cx), np.where(gone, np.nan, cy)   # NaN: never in conflict
    dls = dlat_static(cx, cy, h, v, static_obs) if method == 'mmts' else None
    return Traj(cx, cy, h, v, dls, 'veh', s=s, vprof=v, delta=mode['delta'])


def ncv_hypotheses(nc, xnow, vel, t, static_obs, method):
    """Predicted trajectories of a non-connected vehicle under the three hypotheses."""
    n = int(round(P['T'] / P['dt_chk'])) + 1
    tau = np.arange(n) * P['dt_chk']
    c, s = math.cos(nc.lane_psi), math.sin(nc.lane_psi)
    v_along = vel[0] * c + vel[1] * s
    v_lat = -vel[0] * s + vel[1] * c
    cx0, cy0 = xnow[0] + OFF * math.cos(xnow[2]), xnow[1] + OFF * math.sin(xnow[2])
    out = []
    for k in range(3):
        if k == 0:     # cautious: decelerate, keep lane
            va = np.maximum(0, v_along - 1.5 * tau)
            along = np.concatenate([[0], np.cumsum(0.5 * (va[1:] + va[:-1]) * P['dt_chk'])])
            lat = np.zeros(n)
        elif k == 1:   # nominal: constant speed, keep lane
            along, lat, va = v_along * tau, np.zeros(n), np.full(n, v_along)
        else:          # aggressive: continue lateral motion and accelerate
            va = v_along + 1.0 * tau
            along = v_along * tau + 0.5 * tau ** 2
            lat = v_lat * np.minimum(tau, 2.0) if abs(v_lat) > 0.2 else np.zeros(n)
        cx = cx0 + c * along - s * lat
        cy = cy0 + s * along + c * lat
        dpx, dpy = np.gradient(cx), np.gradient(cy)
        psi = np.where(np.hypot(dpx, dpy) > 1e-3, np.arctan2(dpy, dpx), xnow[2])
        far = tau > P['T_unc']
        cx, cy = np.where(far, np.nan, cx), np.where(far, np.nan, cy)
        dls = dlat_static(cx, cy, psi, va, static_obs) if method == 'mmts' else None
        out.append(Traj(cx, cy, psi, va, dls, 'veh', delta=0.0, tag=('ncv', k)))
    return out


def ncv_constvel(xnow, vel, static_obs, method):
    n = int(round(P['T'] / P['dt_chk'])) + 1
    tau = np.arange(n) * P['dt_chk']
    cx0, cy0 = xnow[0] + OFF * math.cos(xnow[2]), xnow[1] + OFF * math.sin(xnow[2])
    cx, cy = cx0 + vel[0] * tau, cy0 + vel[1] * tau
    far = tau > P['T_unc']
    cx, cy = np.where(far, np.nan, cx), np.where(far, np.nan, cy)
    psi = np.full(n, xnow[2])
    v = np.full(n, math.hypot(*vel))
    dls = dlat_static(cx, cy, psi, v, static_obs) if method == 'mmts' else None
    return Traj(cx, cy, psi, v, dls, 'veh', delta=0.0)


def disc_traj(p0, vel, rad):
    n = int(round(P['T'] / P['dt_chk'])) + 1
    tau = np.arange(n) * P['dt_chk']
    return Traj(p0[0] + vel[0] * tau, p0[1] + vel[1] * tau, None, None, None, 'disc', rad=rad, delta=0.0)


# ----------------------------------------------------------------------------
# Pairwise conflict check (Eq. 12)
# ----------------------------------------------------------------------------
def pair_conflict(method, A, B):
    """Earliest sample index (>=1) at which the envelopes of A and B violate h >= 0,
    returned as (index, penetration) or None."""
    delta = max(A.delta, B.delta)
    d = np.hypot(A.cx - B.cx, A.cy - B.cy)
    if method == 'stcs':
        R = P['stcs_R']
        lim = (R + B.rad if B.kind == 'disc' else 2 * R) + P['stcs_delta']
        bad = np.nonzero(d[1:] < lim)[0]
        if len(bad) == 0:
            return None
        m = bad[0] + 1
        return m, lim - d[m]
    if method == 'fecbf':
        lim = (P['r'] + B.rad if B.kind == 'disc' else 2 * P['r']) + delta
        bad = np.nonzero(d[1:] < lim)[0]
        if len(bad) == 0:
            return None
        m = bad[0] + 1
        return m, lim - d[m]
    # M-MTS: morphing superellipse envelopes
    aA = P['r'] * (1 + P['k_par'] * np.maximum(A.v, 0))
    radA = np.hypot(aA, P['r'])
    if B.kind == 'disc':
        radB = B.rad
    else:
        radB = np.hypot(P['r'] * (1 + P['k_par'] * np.maximum(B.v, 0)), P['r'])
    radB = radB[1:] if np.ndim(radB) else radB
    cand = np.nonzero(d[1:] < radA[1:] + radB + delta)[0] + 1
    # exact check in chronological chunks; stop at the first conflicting sample
    for c0 in range(0, len(cand), 8):
        idx = cand[c0:c0 + 8]
        h = h_traj_batch(A, B, idx, delta)
        bad = np.nonzero(h < 0)[0]
        if len(bad):
            return int(idx[bad[0]]), float(-h[bad[0]])
    return None


def superellipse_batch(cx, cy, psi, a, b, n):
    e = (2.0 / n)[:, None]
    x = a[:, None] * np.sign(_ca)[None] * np.abs(_ca)[None] ** e
    y = b[:, None] * np.sign(_sa)[None] * np.abs(_sa)[None] ** e
    c, s = np.cos(psi)[:, None], np.sin(psi)[:, None]
    return np.stack([cx[:, None] + c * x - s * y, cy[:, None] + s * x + c * y], -1)


def poly_sd_batch(PA, PB):
    sep = None
    for A, B in ((PA, PB), (PB, PA)):
        E = np.roll(A, -1, 1) - A
        N = np.stack([E[..., 1], -E[..., 0]], -1)
        N /= np.sqrt((N * N).sum(-1))[..., None]
        proj = np.einsum('bkd,bmd->bkm', B, N).min(1)
        s = (proj - (A * N).sum(-1)).max(1)
        sep = s if sep is None else np.maximum(sep, s)
    out = sep.copy()
    pos = sep > 0
    if np.any(pos):
        best = None
        for A, B in ((PA[pos], PB[pos]), (PB[pos], PA[pos])):
            E = np.roll(B, -1, 1) - B
            d = A[:, :, None, :] - B[:, None, :, :]
            t = np.clip((d * E[:, None]).sum(-1) / (E * E).sum(-1)[:, None], 0, 1)
            pr = d - t[..., None] * E[:, None]
            m = np.sqrt((pr * pr).sum(-1).min((1, 2)))
            best = m if best is None else np.minimum(best, m)
        out[pos] = best
    return out


def poly_disc_sd_batch(PA, c, rad):
    E = np.roll(PA, -1, 1) - PA
    d = c[:, None, :] - PA
    t = np.clip((d * E).sum(-1) / (E * E).sum(-1), 0, 1)
    pr = d - t[..., None] * E
    dist = np.sqrt((pr * pr).sum(-1).min(1))
    inside = np.all(E[..., 0] * d[..., 1] - E[..., 1] * d[..., 0] >= 0, axis=1)
    return np.where(inside, -dist, dist) - rad


def h_traj_batch(A, B, idx, delta):
    """Barrier (8) between the envelopes of A and B at the sample indices idx (vectorised)."""
    ax, ay, ap, av = A.cx[idx], A.cy[idx], A.psi[idx], A.v[idx]
    bx, by = B.cx[idx], B.cy[idx]
    if B.kind == 'disc':
        dl = smin2(A.dls[idx], dlat_partner(ax, ay, ap, av, bx, by, False, B.rad))
        a, b, n = env_params('mmts', av, dl)
        PA = superellipse_batch(ax, ay, ap, a, b, n)
        return poly_disc_sd_batch(PA, np.stack([bx, by], 1), B.rad) - delta
    bp, bv = B.psi[idx], B.v[idx]
    dlA = smin2(A.dls[idx], dlat_partner(ax, ay, ap, av, bx, by))
    dlB = smin2(B.dls[idx], dlat_partner(bx, by, bp, bv, ax, ay))
    a1, b1, n1 = env_params('mmts', av, dlA)
    a2, b2, n2 = env_params('mmts', bv, dlB)
    PA = superellipse_batch(ax, ay, ap, a1, b1, n1)
    PB = superellipse_batch(bx, by, bp, a2, b2, n2)
    return poly_sd_batch(PA, PB) - delta


def h_traj(A, B, m, delta):
    if B.kind == 'disc':
        dl = smin2(A.dls[m], dlat_partner(A.cx[m], A.cy[m], A.psi[m], A.v[m], B.cx[m], B.cy[m], False, B.rad))
        a, b, n = env_params('mmts', A.v[m], dl)
        PA = superellipse(A.cx[m], A.cy[m], A.psi[m], float(a), float(b), float(n))
        return poly_disc_sd(PA, (B.cx[m], B.cy[m]), B.rad) - delta
    dlA = smin2(A.dls[m], dlat_partner(A.cx[m], A.cy[m], A.psi[m], A.v[m], B.cx[m], B.cy[m]))
    dlB = smin2(B.dls[m], dlat_partner(B.cx[m], B.cy[m], B.psi[m], B.v[m], A.cx[m], A.cy[m]))
    a1, b1, n1 = env_params('mmts', A.v[m], dlA)
    a2, b2, n2 = env_params('mmts', B.v[m], dlB)
    PA = superellipse(A.cx[m], A.cy[m], A.psi[m], float(a1), float(b1), float(n1))
    PB = superellipse(B.cx[m], B.cy[m], B.psi[m], float(a2), float(b2), float(n2))
    return poly_sd(PA, PB) - delta


# ----------------------------------------------------------------------------
# DFCS (Algorithm 1)
# ----------------------------------------------------------------------------
class Budget(Exception):
    pass


class Planner:
    def __init__(self, method):
        self.method = method
        self.timing = dict(gen=0.0, search=0.0, select=0.0)
        self.n_cycles = 0
        self.n_fail = 0
        self.n_budget = 0
        self.nodes = []
        self.wait = {}        # cumulative yield time per agent (priority aging)
        self.prev = {}        # previous cycle's configuration (warm start)

    def eff_prio(self, q):
        return self.agents[q].prio + P["aging"] * self.wait.get(q, 0.0)

    # -- candidate cache -----------------------------------------------------
    def cand(self, i, k, d):
        key = (i, k, d)
        c = self._cand.get(key)
        if c is None:
            c = agent_traj(self.agents[i], k, d, self.static_obs, self.method)
            self._cand[key] = c
        return c

    def pair(self, ki, kj):
        key = (ki, kj) if str(ki) < str(kj) else (kj, ki)
        r = self._pair.get(key, 0)
        if r == 0:
            A = self.cand(*ki) if ki[0] != 'u' else self.unc[ki[1]]
            B = self.cand(*kj) if kj[0] != 'u' else self.unc[kj[1]]
            if A.kind == 'disc':
                A, B = B, A
            r = pair_conflict(self.method, A, B)
            self._pair[key] = r
        return r

    def earliest(self, cfg):
        best = None
        act = self.active
        for a_ in range(len(act)):
            i = act[a_]
            ki = (i,) + cfg[i]
            for b_ in range(a_ + 1, len(act)):
                j = act[b_]
                r = self.pair(ki, (j,) + cfg[j])
                if r is not None and (best is None or r[0] < best[0] or (r[0] == best[0] and r[1] > best[1])):
                    best = (r[0], r[1], i, j)
            for u in range(len(self.unc)):
                r = self.pair(ki, ('u', u))
                if r is not None and (best is None or r[0] < best[0] or (r[0] == best[0] and r[1] > best[1])):
                    best = (r[0], r[1], i, None)
                    self._cu = u
        return best

    def options(self, c, cfg):
        _, _, i, j = c
        agents = [i] if j is None else sorted([i, j], key=self.eff_prio)  # low (effective) priority first
        nd = int(round(P['Dmax'] / P['Delta']))
        acts = []
        for q in agents:
            k, d = cfg[q]
            for d2 in range(d + 1, nd + 1):
                acts.append((q, (k, d2)))
            if self.method != 'stcs':
                for k2 in (0, 1, 2):
                    if k2 != k:
                        acts.append((q, (k2, d)))
        # all actions are kept (completeness w.r.t. the action set); those that remove the
        # current conflict are tried first, in the original order (yield before switch,
        # lower-priority agent before higher-priority agent)
        other = {i: j, j: i}

        def resolves(a):
            q, st = a
            o = other[q]
            ko = ('u', self._cu) if o is None else (o,) + cfg[o]
            return self.pair((q,) + st, ko) is None
        good = [a for a in acts if resolves(a)]
        rest = [a for a in acts if not resolves(a)]
        return good + rest

    def dfcs(self, cfg, visited):
        key = tuple(cfg[i] for i in self.active)
        if key in visited:
            return None
        if self._nodes >= P['node_cap']:
            raise Budget
        visited.add(key)
        self._nodes += 1
        c = self.earliest(cfg)
        if c is None:
            return cfg
        for q, st in self.options(c, cfg):
            cfg2 = dict(cfg)
            cfg2[q] = st
            R = self.dfcs(cfg2, visited)
            if R is not None:
                return R
        return None

    # -- cost (Section IV-E) ---------------------------------------------------
    def cost(self, i, tr, weights):
        C = P['cost']
        ag = self.agents[i]
        rem = ag.path.S - tr.s[-1]
        vmode = max(0.5, float(tr.vprof.max()))
        acc = np.diff(tr.vprof) / P['dt_chk']
        jerk = np.diff(acc) / P['dt_chk']
        J0 = (C['goal'] * rem + C['time'] * (P['T'] + rem / vmode) + C['effort'] * (acc ** 2).sum() * P['dt_chk']
              + C['jerk'] * (jerk ** 2).sum() * P['dt_chk'])
        Jk = np.zeros(3)
        for u, U in enumerate(self.unc):
            d = np.nan_to_num(np.hypot(tr.cx - U.cx, tr.cy - U.cy), nan=1e6) - 2 * P['r']
            pen = C['clear'] * (np.maximum(0, C['c_ref'] - d) ** 2).sum() * P['dt_chk']
            if U.tag is not None and U.tag[0] == 'ncv':
                Jk[U.tag[1]] += pen
            else:
                Jk += pen
        return J0 + Jk

    # -- one planning cycle ------------------------------------------------------
    def plan(self, agents, static_obs, unc, unc_weights):
        """agents: list of Agent; unc: list of Traj (non-modifiable); unc_weights: per-hypothesis weights."""
        self.agents, self.static_obs, self.unc = agents, static_obs, unc
        self.unc_w = unc_weights
        self._cand, self._pair = {}, {}
        self.active = [a.i for a in agents if not a.done]
        modes = [NOMINAL] if self.method == 'stcs' else [0, 1, 2]
        t0 = time.perf_counter()
        # (i) candidate generation and envelope construction, initial preference
        cfg = {}
        for i in self.active:
            best = None
            for k in modes:
                tr = self.cand(i, k, 0)
                J = float(self.cost(i, tr, None) @ unc_weights) if self.method != 'stcs' else 0.0
                if best is None or J < best[0]:
                    best = (J, k)
            cfg[i] = (best[1], 0)
        # plan reuse: last cycle's configuration with its delays shifted by one period
        shift = int(round(P['plan_period'] / P['Delta']))
        warm = {i: (self.prev[i][0], max(0, self.prev[i][1] - shift)) for i in self.active if i in self.prev}
        t1 = time.perf_counter()
        # (ii) depth-first conflict search
        self._nodes = 0
        fail = False
        if len(warm) == len(self.active) and self.earliest(warm) is None:
            R = warm
        else:
            try:
                R = self.dfcs(dict(cfg), set())
                if R is None:
                    fail = True
            except Budget:
                R = None
                fail = True
                self.n_budget += 1
        if fail:
            # fallback: every agent that is (transitively) in conflict brakes and stops
            self.n_fail += 1
            stop = (0 if self.method != 'stcs' else NOMINAL, int(round(P['Dmax'] / P['Delta'])))
            nd = int(round(P['Dmax'] / P['Delta']))
            # greedy prioritised assignment: each agent, in order of decreasing priority, takes the
            # first (smallest-delay) candidate that is conflict-free w.r.t. the agents already assigned
            R = {}
            for i in sorted(self.active, key=lambda q: -self.eff_prio(q)):
                choice = stop
                for d in range(nd + 1):
                    for k in ([cfg[i][0]] + [m for m in modes if m != cfg[i][0]]):
                        ok = all(self.pair((i, k, d), (j,) + R[j]) is None for j in R) and \
                             all(self.pair((i, k, d), ('u', u)) is None for u in range(len(self.unc)))
                        if ok:
                            choice = (k, d)
                            break
                    if choice != stop:
                        break
                R[i] = choice
            changed = True
            while changed:
                changed = False
                for a_ in range(len(self.active)):
                    i = self.active[a_]
                    for j in self.active[a_ + 1:]:
                        if self.pair((i,) + R[i], (j,) + R[j]) is not None:
                            for q in (i, j):
                                if R[q] != stop:
                                    R[q] = stop
                                    changed = True
                    for u in range(len(self.unc)):
                        if R[i] != stop and self.pair((i,) + R[i], ('u', u)) is not None:
                            R[i] = stop
                            changed = True
        self.nodes.append(self._nodes)
        t2 = time.perf_counter()
        # (iii) selection among conflict-free candidates (Eq. 14)
        if self.method != 'stcs' and not fail:
            for i in sorted(self.active, key=lambda q: -self.eff_prio(q)):
                k0, d0 = R[i]
                best = None
                for k in modes:
                    ok = True
                    ki = (i, k, d0)
                    for j in self.active:
                        if j != i and self.pair(ki, (j,) + R[j]) is not None:
                            ok = False
                            break
                    if ok:
                        for u in range(len(self.unc)):
                            if self.pair(ki, ('u', u)) is not None:
                                ok = False
                                break
                    if not ok:
                        continue
                    J = float(self.cost(i, self.cand(i, k, d0), None) @ unc_weights)
                    if best is None or J < best[0]:
                        best = (J, k)
                if best is not None:
                    R[i] = (best[1], d0)
        t3 = time.perf_counter()
        self.timing['gen'] += t1 - t0
        self.timing['search'] += t2 - t1
        self.timing['select'] += t3 - t2
        self.n_cycles += 1
        out = {}
        for i in self.active:
            k, d = R[i]
            if d > 0:
                self.wait[i] = self.wait.get(i, 0.0) + P['plan_period']
            self.prev[i] = (k, d)
            tr = self.cand(i, k, d)
            out[i] = dict(s=tr.s, v=tr.vprof, mode=k)
        return out, t3 - t0


# ----------------------------------------------------------------------------
# CBF-QP (Eqs. 10-11)
# ----------------------------------------------------------------------------
def f_vec(x):
    return np.array([x[3] * math.cos(x[2]), x[3] * math.sin(x[2]), 0.0, 0.0])


def env_poly_state(method, x, ctx_static, partner=None):
    cx, cy = centre(x)
    psi, v = x[2], x[3]
    if method == 'fecbf':
        return superellipse(cx, cy, psi, P['r'], P['r'], 2.0)
    dl = dlat_static(np.array([cx]), np.array([cy]), np.array([psi]), np.array([v]), ctx_static)[0]
    if partner is not None:
        px, py, shared, orad = partner
        dl = float(smin2(dl, dlat_partner(cx, cy, psi, v, px, py, shared, orad)))
    a, b, n = env_params('mmts', v, dl)
    return superellipse(cx, cy, psi, float(a), float(b), float(n))


ALPHA = 4.0   # braking capability assumed by the FE-CBF barrier certificate [m/s^2]


def h_wang(xi, cj, vj, Ds, alpha):
    """Safety barrier certificate of Wang et al. for circular agents with bounded braking."""
    ci = centre(xi)
    vi = (xi[3] * math.cos(xi[2]), xi[3] * math.sin(xi[2]))
    dx, dy = ci[0] - cj[0], ci[1] - cj[1]
    dist = max(math.hypot(dx, dy), 1e-6)
    hs = dist - Ds
    root = math.sqrt(2 * alpha * hs) if hs >= 0 else -math.sqrt(2 * alpha * -hs)
    return root + (dx * (vi[0] - vj[0]) + dy * (vi[1] - vj[1])) / dist


def h_pair_state(method, xi, xj, ctx_i, ctx_j, delta, recip=True):
    if method == 'fecbf':
        cj = centre(xj)
        vj = (xj[3] * math.cos(xj[2]), xj[3] * math.sin(xj[2]))
        return h_wang(xi, cj, vj, 2 * P['r'] + delta, 2 * ALPHA if recip else ALPHA)
    cj = centre(xj)
    ci = centre(xi)
    Pi = env_poly_state(method, xi, ctx_i, (cj[0], cj[1], True, 0.0))
    Pj = env_poly_state(method, xj, ctx_j, (ci[0], ci[1], True, 0.0))
    return poly_sd(Pi, Pj) - delta


def h_disc_state(method, xi, ctx_i, c, rad, delta, vel=None):
    if method == 'fecbf':
        return h_wang(xi, c, (0.0, 0.0) if vel is None else vel, P['r'] + rad + delta, ALPHA)
    Pi = env_poly_state(method, xi, ctx_i, (c[0], c[1], False, rad))
    return poly_disc_sd(Pi, c, rad) - delta


def qp2(uref, A, b, lo, hi):
    """min (u-uref)^T W (u-uref) s.t. A u >= b, lo <= u <= hi, W = diag(1, w_kappa). Exact 2-D active-set enumeration."""
    sc = np.array([1.0, math.sqrt(P['w_kappa'])])
    zr = uref * sc
    rows = [r / sc for r in A] + [np.array([1.0, 0]), np.array([-1.0, 0]), np.array([0, 1.0]), np.array([0, -1.0])]
    rhs = list(b) + [lo[0] * sc[0], -hi[0] * sc[0], lo[1] * sc[1], -hi[1] * sc[1]]
    G = np.array(rows)
    h = np.array(rhs)
    nrm = np.sqrt((G * G).sum(1))
    keep = nrm > 1e-9
    # rows that u cannot influence are dropped; if one of them is violated the QP is infeasible
    stuck = bool(np.any(~keep & (h > 1e-9)))
    G, h = G[keep], h[keep]
    feas = lambda z: np.all(G @ z >= h - 1e-7)
    if feas(zr):
        return zr / sc, not stuck
    cands = []
    viol = G @ zr - h
    for i in range(len(G)):
        if viol[i] < 0:
            cands.append(zr - viol[i] / (G[i] @ G[i]) * G[i])
    m = len(G)
    for i in range(m):
        for j in range(i + 1, m):
            Mx = np.array([G[i], G[j]])
            det = Mx[0, 0] * Mx[1, 1] - Mx[0, 1] * Mx[1, 0]
            if abs(det) > 1e-12:
                cands.append(np.linalg.solve(Mx, np.array([h[i], h[j]])))
    best = None
    for z in cands:
        if feas(z):
            dd = (z - zr) @ (z - zr)
            if best is None or dd < best[0]:
                best = (dd, z)
    if best is None:
        return None, False
    return best[1] / sc, not stuck


def cbf_filter(method, i, x, uref, delta, neighbours, ctx_i, t):
    """neighbours: list of dict(kind, x/c, vel, ctx, delta, recip)."""
    gamma = P['gamma']
    A, b = [], []
    eps = 1e-3
    fx = f_vec(x)
    for nb in neighbours:
        if nb['kind'] in ('veh',):
            xj = nb['x']
            hf = lambda xx, xj=xj: h_pair_state(method, xx, xj, ctx_i, nb['ctx'], max(delta, nb['delta']), nb['recip'])
            h0 = hf(x)
            dh_f = (hf(x + eps * fx) - h0) / eps
            dh_psi = (hf(x + np.array([0, 0, eps, 0])) - h0) / eps
            dh_v = (hf(x + np.array([0, 0, 0, eps])) - h0) / eps
            if nb['recip']:
                g = 0.5 * gamma
                tterm = 0.0
            else:
                g = gamma
                vj = nb['vel']
                xj2 = xj.copy()
                xj2[0] += eps * vj[0]
                xj2[1] += eps * vj[1]
                tterm = (h_pair_state(method, x, xj2, ctx_i, nb['ctx'], max(delta, nb['delta']), False) - h0) / eps
        else:  # disc obstacle (static or dynamic)
            c, rad, vel = nb['c'], nb['rad'], nb['vel']
            hf = lambda xx, c=c: h_disc_state(method, xx, ctx_i, c, rad, delta, vel)
            h0 = hf(x)
            dh_f = (hf(x + eps * fx) - h0) / eps
            dh_psi = (hf(x + np.array([0, 0, eps, 0])) - h0) / eps
            dh_v = (hf(x + np.array([0, 0, 0, eps])) - h0) / eps
            g = gamma
            if vel is not None and (vel[0] != 0 or vel[1] != 0):
                c2 = (c[0] + eps * vel[0], c[1] + eps * vel[1])
                tterm = (h_disc_state(method, x, ctx_i, c2, rad, delta, vel) - h0) / eps
            else:
                tterm = 0.0
        # dh/dx (f + g u) + tterm >= -g h     with  g(x) u = [0, 0, v kappa, a]
        A.append(np.array([dh_v, dh_psi * x[3]]))
        b.append(-g * h0 - dh_f - tterm)
    lo = np.array([max(P['a_min'], -x[3] / P['dt']), -P['kappa_max']])
    hi = np.array([P['a_max'], P['kappa_max']])
    if not A:
        return np.clip(uref, lo, hi), True
    u, ok = qp2(np.clip(uref, lo, hi), A, b, lo, hi)
    if u is None:   # no input satisfies all constraints: brake
        return np.array([lo[0], float(np.clip(uref[1], lo[1], hi[1]))]), False
    return u, ok
