"""Fig. 3 (scenario layouts) and Fig. 7 (N1 interaction) drawn from simulation runs."""
import sys, os, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim, mmts

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.size': 8})


def draw_vehicle(ax, x, env=None, color='#2f6db3', alpha=0.85):
    c = mmts.centre(x)
    if env is not None:
        ax.fill(env[:, 0], env[:, 1], color=color, alpha=0.18, lw=0)
        ax.plot(np.r_[env[:, 0], env[:1, 0]], np.r_[env[:, 1], env[:1, 1]], color=color, lw=0.6)
    F = mmts.rect(c[0], c[1], x[2])
    ax.fill(F[:, 0], F[:, 1], color=color, alpha=alpha, lw=0)


def fig3(seed=1000):
    fig, axes = plt.subplots(3, 4, figsize=(3.5, 3.3))
    fam = ['F', 'S', 'D', 'N']
    titles = ['(a) Free', '(b) Static', '(c) Dynamic', '(d) Non-conn.']
    for col, f in enumerate(fam):
        for row in range(3):
            name = '%s%d' % (f, row + 1)
            ax = axes[row, col]
            sc = sim.make_scenario(name, seed)
            for o in sc['obs']:
                ax.add_patch(plt.Circle(o[:2], o[2], color='#c0392b', lw=0))
            for p in sc['paths']:
                ax.plot(p.xy[:, 0], p.xy[:, 1], color='#2f6db3', lw=0.5, ls='--')
            for i, (s, g, v0) in enumerate(sc['agents']):
                p = sc['paths'][i]
                x = np.array([p.xy[0, 0], p.xy[0, 1], p.h[0], v0])
                c = mmts.centre(x)
                near = sc['obs'][np.hypot(sc['obs'][:, 0] - c[0], sc['obs'][:, 1] - c[1]) < 15] if len(sc['obs']) else sc['obs']
                draw_vehicle(ax, x, mmts.env_poly_state('mmts', x, near))
                ax.plot(*g, marker='x', color='#2f6db3', ms=2.5, mew=0.6)
            for d in sc['dyn']:
                p0 = d.pos(0.0)
                ax.add_patch(plt.Circle(p0, d.rad, color='#c0392b', alpha=0.55, lw=0))
                v = d.vel / max(np.linalg.norm(d.vel), 1e-9) * 6
                ax.arrow(p0[0], p0[1], v[0], v[1], head_width=1.4, color='#c0392b', lw=0.8, length_includes_head=True)
            for fn, lp in sc['ncv']:
                xs = np.array(fn(0.0))
                draw_vehicle(ax, xs, None, color='#7f7f7f', alpha=0.9)
                tr = np.array([mmts.centre(np.array(fn(t))) for t in np.linspace(0, 8, 40)])
                ax.plot(tr[:, 0], tr[:, 1], color='#7f7f7f', lw=0.8)
                ax.annotate('', xy=tr[-1], xytext=tr[-4], arrowprops=dict(arrowstyle='->', color='#7f7f7f', lw=0.8))
            allp = np.concatenate([p.xy for p in sc['paths']])
            cx, cy = allp[:, 0].mean(), allp[:, 1].mean()
            R = max(np.ptp(allp[:, 0]), np.ptp(allp[:, 1])) / 2 + 4
            ax.set_xlim(cx - R, cx + R)
            ax.set_ylim(cy - R, cy + R)
            ax.set_aspect('equal')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(name, fontsize=7, pad=2)
            for sp in ax.spines.values():
                sp.set_linewidth(0.5)
        axes[2, col].set_xlabel(titles[col], fontsize=7)
    fig.tight_layout(h_pad=0.3, w_pad=0.2, pad=0.2)
    fig.savefig(os.path.join(FIG, 'fig3_scenarios.png'), dpi=300)
    plt.close(fig)


def fig7(seed=1000):
    sc = sim.make_scenario('N1', seed)
    r = sim.run_trial('mmts', sc, seed, trace=True)
    tr = r['trace']
    t = np.array(tr['t'])
    ncv = np.array([s[0] for s in tr['ncv']])
    W = np.array([w[0] for w in tr['w']])
    X = np.array(tr['X'])
    lat = ncv[:, 1]
    start = np.argmax(np.abs(lat - lat[0]) > 0.05)
    tc = t[start]
    snaps = [max(0, start - 8), start + 9, start + 45]
    fig = plt.figure(figsize=(3.5, 3.9))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 1.0, 1.0])
    for j, k in enumerate(snaps):
        ax = fig.add_subplot(gs[0, j])
        for y in (-3.5, 3.5):
            ax.axhline(y, color='k', lw=0.6)
        ax.axhline(0, color='k', lw=0.3, ls='--')
        for i in range(X.shape[1]):
            draw_vehicle(ax, X[k, i], tr['env'][k][i])
        draw_vehicle(ax, ncv[k], None, color='#7f7f7f')
        c_ = mmts.centre(ncv[k])
        if k > 5:
            prev = mmts.centre(ncv[k - 5])
            vel = (np.array(c_) - np.array(prev)) / (t[k] - t[k - 5])
            for h, colr in zip(range(3), ['#3a9a5b', '#7f7f7f', '#c0392b']):
                pts = []
                for tau in np.linspace(0, 3, 13):
                    al, la = sim.hyp_offset(h, vel[0], vel[1], tau)
                    pts.append((c_[0] + al, c_[1] + la))
                pts = np.array(pts)
                ax.plot(pts[:, 0], pts[:, 1], color=colr, lw=0.3 + 2.0 * W[k, h], alpha=0.9)
        ego = mmts.centre(X[k, 0])
        ax.set_xlim(ego[0] - 6, ego[0] + 18)
        ax.set_ylim(-5, 5)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('(%s) %.1f s' % ('abc'[j], t[k]), fontsize=7, pad=2)
    ax = fig.add_subplot(gs[1, :])
    for h, lab, colr in zip(range(3), ['cautious', 'nominal', 'aggressive'], ['#3a9a5b', '#7f7f7f', '#c0392b']):
        ax.plot(t, W[:, h], color=colr, lw=1.1, label=lab)
    ax.axvline(tc, color='k', ls=':', lw=0.8)
    ax.set_ylabel('Mode weight', fontsize=7)
    ax.set_xlim(0, min(t[-1], tc + 10))
    ax.tick_params(labelsize=6.5)
    ax.legend(frameon=False, ncol=1, fontsize=6.5, loc='center right')
    ax.text(0.01, 0.92, '(d)', transform=ax.transAxes, fontsize=7, va='top')
    ax = fig.add_subplot(gs[2, :])
    ax.plot(t, X[:, 0, 3], color='#2f6db3', lw=1.1, label='ego speed (m/s)')
    gap = [mmts.poly_sd(mmts.rect(*mmts.centre(X[k, 0]), X[k, 0, 2]), mmts.rect(*mmts.centre(ncv[k]), ncv[k][2])) for k in range(len(t))]
    ax.plot(t, gap, color='#c0392b', lw=1.0, label='gap to NCV (m)')
    ax.axvline(tc, color='k', ls=':', lw=0.8)
    ax.set_xlim(0, min(t[-1], tc + 10))
    ax.set_ylim(0, 12)
    ax.set_xlabel('Time (s)', fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    ax.text(0.01, 0.92, '(e)', transform=ax.transAxes, fontsize=7, va='top')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig7_n1.png'), dpi=300)
    plt.close(fig)
    return dict(t_cutin=float(tc), min_gap=float(np.min(gap)), success=r['success'],
                w_aggr_peak=float(W[:, 2].max()), t_w_aggr_05=float(t[np.argmax(W[:, 2] > 0.5)]) if (W[:, 2] > 0.5).any() else None)


if __name__ == '__main__':
    fig3()
    print(fig7())
