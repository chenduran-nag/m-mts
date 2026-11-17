"""Scenarios and closed-loop simulation for the M-MTS experiments."""
import math, time
import numpy as np
import mmts
from mmts import P, MODES, NOMINAL, Agent, NCV, DynObs, Planner, Traj

LANE = 1.75


# ----------------------------------------------------------------------------
# Scenario helpers
# ----------------------------------------------------------------------------
def right_shift(s, g, off=LANE):
    s, g = np.array(s, float), np.array(g, float)
    d = (g - s) / np.linalg.norm(g - s)
    n = np.array([d[1], -d[0]])
    return tuple(s + off * n), tuple(g + off * n)


def wall_y(x, y0, y1, gaps, r=0.5, step=0.9):
    """Vertical wall of discs at x from y0 to y1 with gaps [(centre, width)] (surface to surface)."""
    out = []
    for y in np.arange(y0, y1 + 1e-6, step):
        if all(abs(y - gc) >= gw / 2 + r for gc, gw in gaps):
            out.append((x, y, r))
    for gc, gw in gaps:
        out.append((x, gc - gw / 2 - r, r))
        out.append((x, gc + gw / 2 + r, r))
    return out


def wall_x(y, x0, x1, r=0.5, step=0.9):
    return [(x, y, r) for x in np.arange(x0, x1 + 1e-6, step)]


def lane_change_fn(x0, y0, v, psi, t_c, dy, dur, v_after=None, decel=1.0):
    """Scripted vehicle moving along heading psi; at t_c it shifts laterally by dy over dur seconds."""
    c, s = math.cos(psi), math.sin(psi)

    def fn(t):
        if v_after is None or t < t_c + dur:
            along = v * t
            vv = v
        else:
            t1 = t - (t_c + dur)
            tb = max(0.0, (v - v_after) / decel)
            if t1 < tb:
                along = v * (t_c + dur) + v * t1 - 0.5 * decel * t1 ** 2
                vv = v - decel * t1
            else:
                along = v * (t_c + dur) + v * tb - 0.5 * decel * tb ** 2 + v_after * (t1 - tb)
                vv = v_after
        if t < t_c:
            lat, vlat = 0.0, 0.0
        elif t < t_c + dur:
            u = (t - t_c) / dur
            lat = dy * (3 * u ** 2 - 2 * u ** 3)
            vlat = dy * (6 * u - 6 * u ** 2) / dur
        else:
            lat, vlat = dy, 0.0
        px = x0 + c * along - s * lat
        py = y0 + s * along + c * lat
        hd = psi + math.atan2(vlat, max(vv, 0.1))
        # convert footprint-centre trajectory to rear-axle state
        return (px - mmts.OFF * math.cos(hd), py - mmts.OFF * math.sin(hd), hd, math.hypot(vv, vlat))
    return fn


def weave_fn(x0, y0, v, psi, amp, period, t0):
    c, s = math.cos(psi), math.sin(psi)
    w = 2 * math.pi / period

    def fn(t):
        tt = max(0.0, t - t0)
        lat = -amp * (1 - math.cos(w * tt)) / 2
        vlat = -amp * w * math.sin(w * tt) / 2 if t > t0 else 0.0
        along = v * t
        px = x0 + c * along - s * lat
        py = y0 + s * along + c * lat
        hd = psi + math.atan2(vlat, v)
        return (px - mmts.OFF * math.cos(hd), py - mmts.OFF * math.sin(hd), hd, math.hypot(v, vlat))
    return fn


def straight_fn(x0, y0, v, psi, t0=0.0):
    c, s = math.cos(psi), math.sin(psi)

    def fn(t):
        along = v * max(0.0, t - t0)
        px, py = x0 + c * along, y0 + s * along
        return (px - mmts.OFF * c, py - mmts.OFF * s, psi, v if t >= t0 else 0.0)
    return fn


# ----------------------------------------------------------------------------
# The twelve scenarios (Fig. 3) and the scaling world
# ----------------------------------------------------------------------------
SCENARIOS = ['F1', 'F2', 'F3', 'S1', 'S2', 'S3', 'D1', 'D2', 'D3', 'N1', 'N2', 'N3']
NAMES = {
    'F1': 'Opposing lanes', 'F2': 'Perpendicular streams', 'F3': 'Circle swap',
    'S1': 'Narrow Gate', 'S2': 'Obstacle field', 'S3': 'Walled corridor',
    'D1': 'Crossing obstacles', 'D2': 'Junction with obstacle', 'D3': 'Obstacle stream',
    'N1': 'Aggressive Cut-off', 'N2': 'Red-light runner', 'N3': 'Weaving vehicle',
}


def build(name, rng, N=None):
    U = rng.uniform
    sc = dict(agents=[], obs=[], dyn=[], ncv=[], bounds=(-45, 45, -45, 45), T_lim=40.0)
    A = sc['agents']
    if name == 'F1':
        x = U(-2, 2)
        A.append(right_shift((-25 + x, 0), (25, 0)) + (0.0,))
        A.append(right_shift((25 + U(-2, 2), 0), (-25, 0)) + (0.0,))
        A.append(right_shift((-35 + x, 0), (15, 0)) + (0.0,))
    elif name == 'F2':
        g = U(8, 12)
        A.append(right_shift((-25, 0), (25, 0)) + (0.0,))
        A.append(right_shift((-25 - g, 0), (25 - g, 0)) + (0.0,))
        o = U(-4, 4)
        A.append(right_shift((0, -25 + o), (0, 25)) + (0.0,))
        A.append(right_shift((0, -25 - g + o), (0, 25 - g)) + (0.0,))
    elif name == 'F3':
        rot = U(0, math.pi / 2)
        for k in range(4):
            th = rot + k * math.pi / 2
            r0 = 20 + U(-3, 3)
            A.append(right_shift((r0 * math.cos(th), r0 * math.sin(th)), (-20 * math.cos(th), -20 * math.sin(th))) + (0.0,))
    elif name == 'S1':
        gc, gw = U(-2, 2), U(3.6, 4.0)
        sc['obs'] = wall_y(0.0, -20, 20, [(gc, gw)])
        ys = [-8 + U(-0.5, 0.5), U(-0.5, 0.5), 8 + U(-0.5, 0.5)]
        for y in ys:
            A.append(((-18.0, y), (18.0, y + U(-2, 2)), 0.0))
        sc['gap'] = (gc, gw)
    elif name == 'S2':
        ys = [-9, -3, 3, 9]
        for y in ys:
            A.append(((-22.0, y + U(-0.5, 0.5)), (22.0, y + U(-0.5, 0.5)), 0.0))
        obs = []
        while len(obs) < 10:
            o = (U(-9, 9), U(-12, 12), U(0.8, 1.5))
            if all(math.hypot(o[0] - p[0], o[1] - p[1]) > o[2] + p[2] + 3.5 for p in obs):
                obs.append(o)
        sc['obs'] = obs
    elif name == 'S3':
        gw = U(7.0, 7.6)
        sc['obs'] = wall_x(gw / 2 + 0.5, -15, 15) + wall_x(-gw / 2 - 0.5, -15, 15)
        A.append(right_shift((-24 + U(-1, 1), 0), (24, 0)) + (0.0,))
        A.append(right_shift((24 + U(-1, 1), 0), (-24, 0)) + (0.0,))
        A.append(right_shift((-33 + U(-1, 1), 0), (16, 0)) + (0.0,))
        sc['corridor'] = gw
    elif name == 'D1':
        A.append(((-25.0, -1.75), (25.0, -1.75), 0.0))
        A.append(((-25.0 - U(6, 9), -5.25), (25.0, -5.25), 0.0))
        for k in range(2):
            x0 = U(-6, 6)
            v = U(2.0, 4.0)
            t_arr = U(2.5, 6.0)              # time at which the obstacle crosses the lanes
            sc['dyn'].append(DynObs(x0, -3.5 - v * t_arr * (1 if k == 0 else -1), 0, v * (1 if k == 0 else -1), 1.0))
    elif name == 'D2':
        A.append(right_shift((-25, 0), (25, 0)) + (0.0,))
        A.append(right_shift((0, -25), (0, 25)) + (0.0,))
        A.append(right_shift((-25 - U(8, 12), 0), (15, 0)) + (0.0,))
        v = U(2.0, 3.5)
        t_arr = U(3.0, 6.0)
        d = v * t_arr / math.sqrt(2)
        sc['dyn'].append(DynObs(d, -d, -v / math.sqrt(2), v / math.sqrt(2), 1.2))
    elif name == 'D3':
        for y in (-6, 0, 6):
            A.append(((-25.0, y + U(-0.5, 0.5)), (25.0, y + U(-0.5, 0.5)), 0.0))
        for k in range(3):
            x0 = U(-10, 10)
            v = U(1.5, 3.5) * (1 if k % 2 == 0 else -1)
            t_arr = U(1.5, 7.0)
            sc['dyn'].append(DynObs(x0, -v * t_arr, 0, v, 1.0))
    elif name == 'N1':
        A.append(((-30.0, -1.75), (60.0, -1.75), 6.0))
        A.append(((-42.0, -1.75), (48.0, -1.75), 6.0))
        x0 = -24 + U(0, 4)
        sc['ncv'].append((lane_change_fn(x0, 1.75, 5.0, 0.0, U(1.0, 2.5), -3.5, 1.5, v_after=3.0, decel=1.0), 0.0))
        sc['bounds'] = (-50, 70, -20, 20)
        sc['T_lim'] = 50.0
    elif name == 'N2':
        A.append(((-25.0, -1.75), (25.0, -1.75), 0.0))
        A.append(((-25.0 - U(8, 12), -1.75), (15.0, -1.75), 0.0))
        v = U(5.0, 7.0)
        t_arr = U(3.5, 7.0)
        x0 = U(-2, 2)
        sc['ncv'].append((straight_fn(x0, -1.75 - v * t_arr, v, math.pi / 2), math.pi / 2))
    elif name == 'N3':
        A.append(((-25.0, -1.75), (35.0, -1.75), 0.0))
        A.append(((-25.0 - U(6, 9), -5.25), (35.0, -5.25), 0.0))
        sc['ncv'].append((weave_fn(-18 + U(0, 6), 1.75, 4.0, 0.0, U(2.5, 3.5), U(5.0, 8.0), U(0.5, 2.0)), 0.0))
        sc['T_lim'] = 45.0
    elif name == 'G':   # scaling world: N agents cross a square between opposite edges
        slots = [-22.75 + 6.5 * k for k in range(8)]
        perm = {'W': rng.permutation(8), 'S': rng.permutation(8)}
        cnt = {'W': 0, 'S': 0}
        for k in range(N):
            side = 'W' if k % 2 == 0 else 'S'
            si = perm[side][cnt[side]]
            cnt[side] += 1
            gi = int(np.clip(si + rng.integers(-2, 3), 0, 7))       # goal within two slots
            a_, b_ = slots[si], slots[gi]
            if side == 'W':
                A.append(((-32.0, a_), (32.0, b_), 0.0))
            else:
                A.append(((a_, -32.0), (b_, 32.0), 0.0))
        obs = []
        tries = 0
        while len(obs) < 6 and tries < 1000:
            tries += 1
            o = (U(-15, 15), U(-15, 15), U(0.8, 1.2))
            if all(math.hypot(o[0] - p[0], o[1] - p[1]) > o[2] + p[2] + 5.0 for p in obs):
                obs.append(o)
        sc['obs'] = obs
        sc['bounds'] = (-38, 38, -38, 38)
        sc['T_lim'] = 60.0
    else:
        raise ValueError(name)
    sc['obs'] = np.array(sc['obs'], float).reshape(-1, 3)
    sc['name'] = name
    return sc


def make_scenario(name, seed, N=None):
    """Sample a scenario instance; resample until every agent has a path."""
    rng = np.random.default_rng(seed)
    for _ in range(50):
        sc = build(name, rng, N)
        paths = []
        for (s, g, v0) in sc['agents']:
            p = mmts.astar_path(s, g, sc['obs'], sc['bounds'])
            if p is None:
                break
            paths.append(p)
        else:
            sc['paths'] = paths
            return sc
    raise RuntimeError('no feasible instance')


def free_time(S, v0):
    """Travel time of an agent alone on its path with the fastest (aggressive) profile."""
    mode = MODES[2]
    s, v, t, dt = 0.0, v0, 0.0, 0.05
    while s < S - 1.0 and t < 500:
        vt = min(mode['v'], math.sqrt(max(0.0, 2 * 2.5 * (S - s - 0.3))))
        dv = min(max(vt - v, -4.0 * dt), mode['acc'] * dt)
        v2 = max(0.0, v + dv)
        s += 0.5 * (v + v2) * dt
        v = v2
        t += dt
    return t


# ----------------------------------------------------------------------------
# Closed-loop trial
# ----------------------------------------------------------------------------
def hyp_offset(k, v_along, v_lat, tau):
    if k == 0:
        tb = min(tau, v_along / 1.5) if v_along > 0 else 0.0
        return v_along * tb - 0.75 * tb ** 2, 0.0
    if k == 1:
        return v_along * tau, 0.0
    lat = v_lat * min(tau, 2.0) if abs(v_lat) > 0.2 else 0.0
    return v_along * tau + 0.5 * tau ** 2, lat


def ncv_vel(nc, t):
    a, b = nc.state(t + 0.05), nc.state(max(0.0, t - 0.05))
    ca = mmts.centre(a)
    cb = mmts.centre(b)
    dt = 0.05 + min(t, 0.05)
    return np.array([(ca[0] - cb[0]) / dt, (ca[1] - cb[1]) / dt])


def run_trial(method, sc, seed, trace=False):
    rng = np.random.default_rng(seed + 7919)
    n = len(sc['agents'])
    prio = rng.permutation(n)
    agents = [Agent(i, sc['paths'][i], sc['agents'][i][2], int(prio[i])) for i in range(n)]
    obs = sc['obs']
    ncvs = [NCV(fn, lp) for fn, lp in sc['ncv']]
    dyns = sc['dyn']
    planner = Planner(method)
    t_free = [free_time(a.path.S, a.x[3]) for a in agents]
    dt = P['dt']
    t = 0.0
    next_plan = 0.0
    min_clear = 1e9
    collision = False
    qp_steps = qp_bad = 0
    runtimes = []
    qp_times = []
    tr = dict(t=[], X=[], ncv=[], env=[], mode=[], w=[], dyn=[], plans=[]) if trace else None
    nplan = int(round(P['T'] / P['dt_chk'])) + 1
    tau_grid = np.arange(nplan) * P['dt_chk']

    def near_obs(c, R=15.0):
        if len(obs) == 0:
            return obs
        d = np.hypot(obs[:, 0] - c[0], obs[:, 1] - c[1])
        return obs[d < R]

    while t < sc['T_lim'] - 1e-9:
        # ---------------- planning ------------------------------------------------
        if t >= next_plan - 1e-9:
            unc, wts = [], np.ones(3) / 3
            for nc in ncvs:
                xs = nc.state(t)
                vel = ncv_vel(nc, t)
                if method == 'stcs':
                    unc.append(mmts.ncv_constvel(xs, vel, obs, method))
                else:
                    hyps = mmts.ncv_hypotheses(nc, xs, vel, t, obs, method)
                    for k in range(3):
                        if nc.w[k] >= P['w_active']:
                            unc.append(hyps[k])
                    wts = nc.w.copy()
            for d in dyns:
                unc.append(mmts.disc_traj(d.pos(t), d.velocity(t), d.rad))
            plans, rt = planner.plan(agents, obs, unc, wts)
            runtimes.append(rt)
            for i, pl in plans.items():
                agents[i].plan = dict(t0=t, s=pl['s'], v=pl['v'])
                agents[i].mode = pl['mode']
            if trace:
                tr['plans'].append((t, {i: (pl['s'].copy(), pl['mode']) for i, pl in plans.items()}))
            next_plan += P['plan_period']
        # ---------------- control ------------------------------------------------------
        U = {}
        ncv_states = [(nc.state(t), ncv_vel(nc, t)) for nc in ncvs]
        for a in agents:
            if a.done:
                continue
            p = np.array([a.x[0], a.x[1]])
            s_now, lat, a.hint = a.path.project(p, a.hint)
            a.s = s_now
            a.track_err.append(lat)
            if s_now >= a.path.S - 1.0:
                a.done, a.t_done = True, t
                a.x[3] = 0.0
                continue
            pl = a.plan
            tau = t - pl['t0']
            s_ref = float(np.interp(tau, tau_grid, pl['s']))
            v_ref = float(np.interp(tau, tau_grid, pl['v']))
            kk = min(int(tau / P['dt_chk']), len(tau_grid) - 2)
            a_ff = (pl['v'][kk + 1] - pl['v'][kk]) / P['dt_chk']      # feed-forward from the profile
            a_ref = a_ff + 1.5 * (v_ref - a.x[3]) + 0.5 * (s_ref - s_now)
            Ld = max(3.0, 0.5 * a.x[3])
            tx, ty, _ = a.path.at(s_now + Ld)
            alpha = math.atan2(float(ty) - a.x[1], float(tx) - a.x[0]) - a.x[2]
            k_ref = 2 * math.sin(alpha) / Ld
            uref = np.array([a_ref, k_ref])
            if method == 'stcs':
                lo = np.array([max(P['a_min'], -a.x[3] / dt), -P['kappa_max']])
                hi = np.array([P['a_max'], P['kappa_max']])
                U[a.i] = np.clip(uref, lo, hi)
                continue
            ci = mmts.centre(a.x)
            ctx_i = near_obs(ci)
            nbs = []
            for b in agents:
                if b is a or b.done:
                    continue
                cb = mmts.centre(b.x)
                if math.hypot(cb[0] - ci[0], cb[1] - ci[1]) > P['sense']:
                    continue
                nbs.append(dict(kind='veh', x=b.x.copy(), ctx=near_obs(cb), delta=MODES[b.mode]['delta'],
                                recip=not b.done, vel=np.zeros(2)))
            for (xs, vel) in ncv_states:
                cb = mmts.centre(xs)
                if math.hypot(cb[0] - ci[0], cb[1] - ci[1]) > P['sense']:
                    continue
                nbs.append(dict(kind='veh', x=xs, ctx=near_obs(cb), delta=0.0, recip=False, vel=vel))
            for o in ctx_i:
                if math.hypot(o[0] - ci[0], o[1] - ci[1]) < 8.0 + o[2]:
                    nbs.append(dict(kind='disc', c=(o[0], o[1]), rad=o[2], vel=None))
            for d in dyns:
                c = d.pos(t)
                if math.hypot(c[0] - ci[0], c[1] - ci[1]) < P['sense']:
                    nbs.append(dict(kind='disc', c=(c[0], c[1]), rad=d.rad, vel=d.velocity(t)))
            q0 = time.perf_counter()
            u, ok = mmts.cbf_filter(method, a.i, a.x, uref, MODES[a.mode]['delta'], nbs, ctx_i, t)
            qp_times.append(time.perf_counter() - q0)
            qp_steps += 1
            qp_bad += (not ok)
            U[a.i] = u
        # ---------------- integrate ----------------------------------------------------
        h = dt / P['substeps']
        for a in agents:
            if a.done or a.i not in U:
                continue
            u = U[a.i]
            for _ in range(P['substeps']):
                x = a.x
                x[0] += h * x[3] * math.cos(x[2])
                x[1] += h * x[3] * math.sin(x[2])
                x[2] += h * x[3] * u[1]
                x[3] = max(0.0, x[3] + h * u[0])
        t += dt
        # ---------------- mode-weight update (Eq. 15) -------------------------------------
        if method != 'stcs' and t >= P['w_window']:
            for nc in ncvs:
                xp = nc.state(t - P['w_window'])
                vp = ncv_vel(nc, t - P['w_window'])
                cpast = mmts.centre(xp)
                cnow = mmts.centre(nc.state(t))
                c_, s_ = math.cos(nc.lane_psi), math.sin(nc.lane_psi)
                va, vl = vp[0] * c_ + vp[1] * s_, -vp[0] * s_ + vp[1] * c_
                lik = np.zeros(3)
                for k in range(3):
                    al, la = hyp_offset(k, va, vl, P['w_window'])
                    px = cpast[0] + c_ * al - s_ * la
                    py = cpast[1] + s_ * al + c_ * la
                    lik[k] = math.exp(-((px - cnow[0]) ** 2 + (py - cnow[1]) ** 2) / (2 * P['sigma'] ** 2))
                w = nc.w * lik
                w = w / w.sum() if w.sum() > 1e-300 else np.ones(3) / 3
                w = np.maximum(w, P['w_floor'])
                nc.w = w / w.sum()
        # ---------------- safety bookkeeping (physical footprints) -----------------------
        foot = [mmts.rect(*mmts.centre(a.x), a.x[2]) for a in agents]
        cents = [mmts.centre(a.x) for a in agents]
        coll_with = None
        for i in range(n):
            if agents[i].done:
                continue
            for j in range(i + 1, n):
                if agents[j].done:
                    continue
                if math.hypot(cents[i][0] - cents[j][0], cents[i][1] - cents[j][1]) < 8.0:
                    d = mmts.poly_sd(foot[i], foot[j])
                    if d <= 0:
                        coll_with = 'agent'
                    min_clear = min(min_clear, d)
                else:
                    min_clear = min(min_clear, math.hypot(cents[i][0] - cents[j][0], cents[i][1] - cents[j][1]) - 2 * mmts.HALF_L)
        for i in range(n):
            if agents[i].done:
                continue
            for o in near_obs(cents[i], 6.0):
                d_ = mmts.poly_disc_sd(foot[i], (o[0], o[1]), o[2])
                coll_with = 'static' if d_ <= 0 else coll_with
                min_clear = min(min_clear, d_)
            for d in dyns:
                c = d.pos(t)
                if math.hypot(c[0] - cents[i][0], c[1] - cents[i][1]) < 6.0 + d.rad:
                    d_ = mmts.poly_disc_sd(foot[i], c, d.rad)
                    coll_with = 'dynamic' if d_ <= 0 else coll_with
                    min_clear = min(min_clear, d_)
            for nc in ncvs:
                xs = nc.state(t)
                cb = mmts.centre(xs)
                if math.hypot(cb[0] - cents[i][0], cb[1] - cents[i][1]) < 8.0:
                    d_ = mmts.poly_sd(foot[i], mmts.rect(cb[0], cb[1], xs[2]))
                    coll_with = 'ncv' if d_ <= 0 else coll_with
                    min_clear = min(min_clear, d_)
        if trace:
            tr['t'].append(t)
            tr['X'].append(np.array([a.x.copy() for a in agents]))
            tr['ncv'].append([nc.state(t) for nc in ncvs])
            tr['w'].append([nc.w.copy() for nc in ncvs])
            tr['mode'].append([a.mode for a in agents])
            tr['dyn'].append([d.pos(t) for d in dyns])
            envs = []
            for a in agents:
                ci = mmts.centre(a.x)
                envs.append(mmts.env_poly_state(method, a.x, near_obs(ci)))
            tr['env'].append(envs)
        if min_clear <= 0:
            collision = True
            break_with = coll_with
            break
        if all(a.done for a in agents):
            break
    success = (not collision) and all(a.done for a in agents)
    stuck = [a.i for a in agents if not a.done]
    res = dict(method=method, scenario=sc['name'], seed=seed, n=n, success=success, collision=collision,
               min_clear=min_clear, runtime_ms=1000 * float(np.mean(runtimes)) if runtimes else float('nan'),
               runtime_max_ms=1000 * float(np.max(runtimes)) if runtimes else float('nan'),
               qp_ms=1000 * float(np.mean(qp_times)) if qp_times else float('nan'),
               qp_infeasible=qp_bad / qp_steps if qp_steps else 0.0,
               dfcs_fail=planner.n_fail / max(1, planner.n_cycles),
               dfcs_budget=planner.n_budget / max(1, planner.n_cycles),
               nodes_mean=float(np.mean(planner.nodes)) if planner.nodes else 0.0,
               timing={k: v / max(1, planner.n_cycles) for k, v in planner.timing.items()},
               t_end=t, coll_with=break_with if collision else None, n_unfinished=len(stuck))
    if success:
        res['subopt'] = sum(a.t_done for a in agents) / sum(t_free)
        res['track_err'] = float(np.mean([np.mean(a.track_err) for a in agents if a.track_err]))
    if trace:
        res['trace'] = tr
    return res
