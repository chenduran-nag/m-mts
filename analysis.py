"""Tables, statistics and figures for the paper, computed from the per-trial results.

Usage: python analysis.py results/scenarios.jsonl results/scaling.jsonl
Writes results/summary.md, results/summary.json and figures/*.png.
"""
import sys, os, json, math, collections
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim, mmts

METHODS = ['mmts', 'fecbf', 'stcs']
LABEL = {'mmts': 'M-MTS', 'fecbf': 'FE-CBF', 'stcs': 'STCS'}
COL = {'mmts': '#2f6db3', 'fecbf': '#d9822b', 'stcs': '#3a9a5b'}
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 8})


def load(p):
    R = [json.loads(l) for l in open(p)]
    bad = [r for r in R if 'error' in r]
    if bad:
        print('WARNING: %d crashed trials' % len(bad), bad[:3])
    return [r for r in R if 'error' not in r]


def ci95(x):
    x = np.asarray(x, float)
    if len(x) < 2:
        return float('nan')
    return float(stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / math.sqrt(len(x)))


def wilson(k, n):
    if n == 0:
        return (float('nan'),) * 2
    z = 1.96
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def fmt(x, c, nd=2):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 'n/a'
    return ('%.' + str(nd) + 'f ± %.' + str(nd) + 'f') % (x, c)


def main():
    RS = load(sys.argv[1])
    RG = load(sys.argv[2]) if len(sys.argv) > 2 and os.path.exists(sys.argv[2]) else []
    by = collections.defaultdict(list)
    for r in RS:
        by[(r['scenario'], r['method'])].append(r)
    out = {}
    md = []

    # ---------- per-scenario table -------------------------------------------------------
    md.append('## Per-scenario results (%d trials per cell)\n' % len(by[('S1', 'mmts')]))
    md.append('| Scenario | Agents | ' + ' | '.join('%s success' % LABEL[m] for m in METHODS) + ' | ' +
              ' | '.join('%s margin (m)' % LABEL[m] for m in METHODS) + ' |')
    md.append('|' + '---|' * (2 + 2 * len(METHODS)))
    per = {}
    for s in sim.SCENARIOS:
        row = [s + ' ' + sim.NAMES[s], str(by[(s, 'mmts')][0]['n'])]
        cells = {}
        for m in METHODS:
            L = by[(s, m)]
            k = sum(r['success'] for r in L)
            cells[m] = dict(k=k, n=len(L), coll=sum(r['collision'] for r in L),
                            margin=float(np.mean([r['min_clear'] for r in L])),
                            margin_ci=ci95([r['min_clear'] for r in L]))
            row.append('%d/%d%s' % (k, len(L), (' (%d coll.)' % cells[m]['coll']) if cells[m]['coll'] else ''))
        for m in METHODS:
            row.append(fmt(cells[m]['margin'], cells[m]['margin_ci']))
        per[s] = cells
        md.append('| ' + ' | '.join(row) + ' |')
    out['per_scenario'] = per

    # ---------- Table I: S1 and N1 in detail ---------------------------------------------
    md.append('\n## Table I: S1 and N1\n')
    tab1 = {}
    for s in ('S1', 'N1'):
        tab1[s] = {}
        for m in METHODS:
            L = by[(s, m)]
            succ = [r for r in L if r['success']]
            e = dict(success_k=sum(r['success'] for r in L), n=len(L),
                     runtime=float(np.mean([r['runtime_ms'] for r in L])), runtime_ci=ci95([r['runtime_ms'] for r in L]),
                     margin=float(np.mean([r['min_clear'] for r in L])), margin_ci=ci95([r['min_clear'] for r in L]),
                     subopt=float(np.mean([r['subopt'] for r in succ])) if succ else float('nan'),
                     subopt_ci=ci95([r['subopt'] for r in succ]) if succ else float('nan'),
                     track=float(np.mean([r['track_err'] for r in succ])) if succ else float('nan'),
                     track_ci=ci95([r['track_err'] for r in succ]) if succ else float('nan'),
                     qp_ms=float(np.nanmean([r['qp_ms'] for r in L])) if m != 'stcs' else float('nan'),
                     qp_inf=float(np.mean([r['qp_infeasible'] for r in L])),
                     dfcs_fail=float(np.mean([r['dfcs_fail'] for r in L])))
            tab1[s][m] = e
        md.append('| Scenario | Metric | ' + ' | '.join(LABEL[m] for m in METHODS) + ' |')
        md.append('|---|---|---|---|---|')
        md.append('| %s | Success | ' % s + ' | '.join('%d/%d' % (tab1[s][m]['success_k'], tab1[s][m]['n']) for m in METHODS) + ' |')
        md.append('| %s | Runtime (ms) | ' % s + ' | '.join(fmt(tab1[s][m]['runtime'], tab1[s][m]['runtime_ci']) for m in METHODS) + ' |')
        md.append('| %s | Safety margin (m) | ' % s + ' | '.join(fmt(tab1[s][m]['margin'], tab1[s][m]['margin_ci']) for m in METHODS) + ' |')
        md.append('| %s | Suboptimality | ' % s + ' | '.join(fmt(tab1[s][m]['subopt'], tab1[s][m]['subopt_ci']) for m in METHODS) + ' |')
        md.append('| %s | Tracking error (m) | ' % s + ' | '.join(fmt(tab1[s][m]['track'], tab1[s][m]['track_ci'], 3) for m in METHODS) + ' |')
    out['table1'] = tab1

    # ---------- statistics -----------------------------------------------------------------
    md.append('\n## Statistical tests\n')
    tests = {}
    for s in sim.SCENARIOS:
        for m in ('fecbf', 'stcs'):
            a, b = by[(s, 'mmts')], by[(s, m)]
            ka, kb = sum(r['success'] for r in a), sum(r['success'] for r in b)
            p_f = stats.fisher_exact([[ka, len(a) - ka], [kb, len(b) - kb]])[1]
            xa, xb = [r['min_clear'] for r in a], [r['min_clear'] for r in b]
            w = stats.ttest_ind(xa, xb, equal_var=False)
            va, vb = np.var(xa, ddof=1) / len(xa), np.var(xb, ddof=1) / len(xb)
            dof = (va + vb) ** 2 / (va ** 2 / (len(xa) - 1) + vb ** 2 / (len(xb) - 1)) if va + vb > 0 else float('nan')
            sa = [r['subopt'] for r in a if r['success']]
            sb = [r['subopt'] for r in b if r['success']]
            if len(sa) > 1 and len(sb) > 1:
                ws = stats.ttest_ind(sa, sb, equal_var=False)
                sub = (float(ws.statistic), float(ws.pvalue))
            else:
                sub = (float('nan'), float('nan'))
            ra, rb = [r['runtime_ms'] for r in a], [r['runtime_ms'] for r in b]
            wr = stats.ttest_ind(ra, rb, equal_var=False)
            tests['%s_%s' % (s, m)] = dict(fisher_p=float(p_f), margin_t=float(w.statistic), margin_p=float(w.pvalue),
                                          margin_df=float(dof), subopt_t=sub[0], subopt_p=sub[1],
                                          runtime_t=float(wr.statistic), runtime_p=float(wr.pvalue))
            md.append('- %s, M-MTS vs %s: success Fisher p = %.3g; margin Welch t = %.2f (df ≈ %.0f), p = %.3g; '
                      'suboptimality t = %.2f, p = %.3g; runtime t = %.2f, p = %.3g' %
                      (s, LABEL[m], p_f, w.statistic, dof, w.pvalue, sub[0], sub[1], wr.statistic, wr.pvalue))
    out['tests'] = tests

    # ---------- clearance distribution (Table / Fig. 8) --------------------------------------
    md.append('\n## Minimum clearance over trials in which the agents drove (success or collision; time-outs excluded)\n')
    clr = {m: np.array([r['min_clear'] for r in RS if r['method'] == m and (r['success'] or r['collision'])]) for m in METHODS}
    pct = {m: {p: float(np.percentile(clr[m], p)) for p in (5, 50, 95)} for m in METHODS}
    for m in METHODS:
        pct[m]['min'] = float(clr[m].min())
        pct[m]['n'] = int(len(clr[m]))
        pct[m]['collisions'] = int(sum(r['collision'] for r in RS if r['method'] == m))
    out['clearance'] = pct
    md.append('| Percentile | ' + ' | '.join(LABEL[m] for m in METHODS) + ' |')
    md.append('|---|---|---|---|')
    for p in (5, 50, 95, 'min'):
        md.append('| %s | ' % p + ' | '.join('%.2f' % pct[m][p] for m in METHODS) + ' |')
    md.append('| collisions | ' + ' | '.join('%d/%d' % (pct[m]['collisions'], pct[m]['n']) for m in METHODS) + ' |')

    # overall
    tot = {}
    for m in METHODS:
        L = [r for r in RS if r['method'] == m]
        tot[m] = dict(success=sum(r['success'] for r in L), n=len(L), collisions=sum(r['collision'] for r in L),
                      timeouts=sum((not r['success']) and (not r['collision']) for r in L),
                      qp_inf=float(np.mean([r['qp_infeasible'] for r in L])),
                      dfcs_fail=float(np.mean([r['dfcs_fail'] for r in L])),
                      dfcs_budget=float(np.mean([r['dfcs_budget'] for r in L])),
                      track=float(np.mean([r['track_err'] for r in L if r['success']])),
                      track_max=float(np.max([r['track_err'] for r in L if r['success']])),
                      runtime=float(np.mean([r['runtime_ms'] for r in L])),
                      runtime_max=float(np.max([r['runtime_max_ms'] for r in L])),
                      qp_ms=float(np.nanmean([r['qp_ms'] for r in L])) if m != 'stcs' else float('nan'))
    out['overall'] = tot
    md.append('\n## Overall\n')
    for m in METHODS:
        md.append('- %s: %s' % (LABEL[m], json.dumps(tot[m])))

    # Fig. 8: KDE of minimum clearance
    fig, ax = plt.subplots(figsize=(3.4, 2.0))
    xs = np.linspace(min(c.min() for c in clr.values()) - 0.3, np.percentile(np.concatenate(list(clr.values())), 99.5) + 0.5, 400)
    for m in METHODS:
        k = stats.gaussian_kde(clr[m])
        ax.fill_between(xs, k(xs), alpha=0.25, color=COL[m])
        ax.plot(xs, k(xs), color=COL[m], lw=1.4, label=LABEL[m])
    ax.axvline(0, color='k', lw=0.6, ls=':')
    ax.set_xlabel('Minimum clearance per trial (m)')
    ax.set_ylabel('Density')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig8_clearance.png'), dpi=300)
    plt.close(fig)

    # ---------- scaling -------------------------------------------------------------------
    if RG:
        sc = collections.defaultdict(list)
        for r in RG:
            sc[(r['N'], r['method'])].append(r)
        Ns = sorted(set(r['N'] for r in RG))
        scal = {}
        md.append('\n## Scaling world\n')
        md.append('| N | ' + ' | '.join('%s success' % LABEL[m] for m in METHODS) + ' | ' + ' | '.join('%s runtime (ms)' % LABEL[m] for m in METHODS) + ' |')
        md.append('|' + '---|' * 7)
        for N in Ns:
            scal[N] = {}
            row = [str(N)]
            for m in METHODS:
                L = sc[(N, m)]
                k = sum(r['success'] for r in L)
                lo, hi = wilson(k, len(L))
                rt = [r['runtime_ms'] for r in L]
                scal[N][m] = dict(k=k, n=len(L), lo=lo, hi=hi, rt=float(np.mean(rt)), rt_ci=ci95(rt),
                                  rt_max=float(np.max([r['runtime_max_ms'] for r in L])),
                                  coll=sum(r['collision'] for r in L),
                                  budget=float(np.mean([r['dfcs_budget'] for r in L])),
                                  fail=float(np.mean([r['dfcs_fail'] for r in L])),
                                  timing={q: float(np.mean([r['timing'][q] for r in L])) for q in ('gen', 'search', 'select')},
                                  qp_ms=float(np.nanmean([r['qp_ms'] for r in L])) if m != 'stcs' else float('nan'))
                row.append('%d/%d' % (k, len(L)) + (' (%d coll.)' % scal[N][m]['coll'] if scal[N][m]['coll'] else ''))
            for m in METHODS:
                row.append(fmt(scal[N][m]['rt'], scal[N][m]['rt_ci'], 1))
            md.append('| ' + ' | '.join(row) + ' |')
        out['scaling'] = {str(k): v for k, v in scal.items()}

        fig, ax = plt.subplots(figsize=(3.4, 2.1))
        for m in METHODS:
            y = np.array([scal[N][m]['rt'] for N in Ns])
            c = np.array([scal[N][m]['rt_ci'] for N in Ns])
            ax.plot(Ns, y, '-o', ms=3, color=COL[m], label=LABEL[m])
            ax.fill_between(Ns, y - c, y + c, color=COL[m], alpha=0.2)
        ax.set_xlabel('Number of agents N')
        ax.set_ylabel('Planning time (ms)')
        ax.set_yscale('log')
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, 'fig4_runtime.png'), dpi=300)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(3.4, 2.1))
        for m in METHODS:
            y = np.array([100 * scal[N][m]['k'] / scal[N][m]['n'] for N in Ns])
            lo = np.array([100 * scal[N][m]['lo'] for N in Ns])
            hi = np.array([100 * scal[N][m]['hi'] for N in Ns])
            ax.plot(Ns, y, '-o', ms=3, color=COL[m], label=LABEL[m])
            ax.fill_between(Ns, lo, hi, color=COL[m], alpha=0.15)
        ax.set_xlabel('Number of agents N')
        ax.set_ylabel('Success rate (%)')
        ax.set_ylim(-3, 103)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, 'fig5_success.png'), dpi=300)
        plt.close(fig)

        # Fig. 6: runtime breakdown of M-MTS
        fig, ax = plt.subplots(figsize=(3.4, 2.0))
        parts = [('gen', 'Candidate generation\n+ envelopes'), ('search', 'Conflict search\n(DFCS)'), ('select', 'Selection')]
        bottom = np.zeros(len(Ns))
        for (q, lab), colr in zip(parts, ['#9ecae1', '#3182bd', '#08519c']):
            y = np.array([1000 * scal[N]['mmts']['timing'][q] for N in Ns])
            ax.bar([str(N) for N in Ns], y, bottom=bottom, color=colr, label=lab.replace('\n', ' '))
            bottom += y
        ax.set_xlabel('Number of agents N')
        ax.set_ylabel('Time per cycle (ms)')
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, 'fig6_breakdown.png'), dpi=300)
        plt.close(fig)

    json.dump(out, open(os.path.join('results', 'summary.json'), 'w'), indent=1, default=float)
    open(os.path.join('results', 'summary.md'), 'w', encoding='utf8').write('\n'.join(md) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(md))


if __name__ == '__main__':
    main()
