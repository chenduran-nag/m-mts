"""Run the experiments of the paper and store per-trial results as JSON lines.

Usage:
  python run_experiments.py scenarios TRIALS OUT.jsonl     # 12 scenarios x 3 methods
  python run_experiments.py scaling   TRIALS OUT.jsonl     # N = 2..16 in the scaling world
"""
import sys, json, os, time
from multiprocessing import Pool
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim

METHODS = ['mmts', 'fecbf', 'stcs']
NS = [2, 4, 6, 8, 10, 12, 14, 16]


def job(args):
    kind, name, N, seed, method = args
    sc = sim.make_scenario(name, seed, N)
    t0 = time.time()
    try:
        r = sim.run_trial(method, sc, seed)
    except Exception as e:  # record crashes instead of hiding them
        return dict(method=method, scenario=name, N=N, seed=seed, error=repr(e))
    r['N'] = N if N is not None else len(sc['agents'])
    r['wall_s'] = time.time() - t0
    return {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in r.items()}


def main():
    kind, trials, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    jobs = []
    if kind == 'scenarios':
        names = sys.argv[4].split(',') if len(sys.argv) > 4 else sim.SCENARIOS
        for name in names:
            for seed in range(trials):
                for m in METHODS:
                    jobs.append((kind, name, None, 1000 + seed, m))
    else:
        for N in NS:
            for seed in range(trials):
                for m in METHODS:
                    jobs.append((kind, 'G', N, 5000 + 100 * N + seed, m))
    # longest jobs first for better load balance
    jobs.sort(key=lambda j: (j[4] != 'mmts', -(j[2] or 0)))
    done = 0
    t0 = time.time()
    with Pool(int(os.environ.get('NPROC', 14))) as pool, open(out, 'w') as f:
        for r in pool.imap_unordered(job, jobs):
            f.write(json.dumps(r) + '\n')
            f.flush()
            done += 1
            if done % 25 == 0 or done == len(jobs):
                print('%d/%d  %.0fs' % (done, len(jobs), time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
