# M-MTS: Adaptive Multi-Modal Trajectory Synthesis with Morphing Safety Envelopes

This repository contains the simulation code for the paper *Adaptive Multi-Modal Trajectory Synthesis with Morphing Safety Envelopes*. It includes:

- the M-MTS planner;
- the two baselines (FE-CBF and STCS);
- the twelve test scenarios and the scaling study;
- the scripts that produce every result table and figure in Section VI of the paper.

The simulation is kinematic and is written in Python.

## Requirements

- Python 3.10 or later
- NumPy, SciPy and Matplotlib (tested with NumPy 2.4, SciPy 1.16 and Matplotlib 3.10)

```
pip install -r requirements.txt
```

## Files

| File | Contents |
|---|---|
| `mmts.py` | The M-MTS method and the code shared with the baselines: envelope law (Eq. 7), signed-distance barrier (Eqs. 8–9), state-lattice path planner, speed profiles, DFCS (Algorithm 1), trajectory selection and the CBF-QP (Eqs. 10–11). All parameters (Table I of the paper) are in `P` and `MODES`. |
| `sim.py` | Scenario definitions (F1–N3 and the scaling world `G`), closed-loop simulation and metrics |
| `run_experiments.py` | Runs all trials in parallel and writes one JSON line per trial |
| `analysis.py` | Statistical tests and result tables (`results/summary.md`, `results/summary.json`), plus Figs. 4–7 |
| `figures.py` | Fig. 3 (scenario layouts) and Fig. 8 (N1 interaction) |
| `results/` | Raw per-trial results used in the paper (`scenarios.jsonl`: 1,800 trials; `scaling.jsonl`: 720 trials) and their summaries |
| `figures/` | The generated figures |

## Reproducing the results

```
python run_experiments.py scenarios 50 results/scenarios.jsonl
python run_experiments.py scaling 30 results/scaling.jsonl
python analysis.py results/scenarios.jsonl results/scaling.jsonl
python figures.py
```

The trial seeds are fixed: 1000–1049 for the scenarios, and 5000 + 100 N + k for the scaling world. Reruns therefore reproduce the same trials. Runtimes depend on the machine. The paper's numbers were measured on an Intel Core i7-10700K with 16 GB RAM, running 14 trials in parallel.

## Mapping to the paper

Some figure file names and the table headings in `results/summary.md` use the numbering from an earlier draft of the paper. They correspond to the paper as follows:

| File / section | Paper |
|---|---|
| `figures/fig3_scenarios.png` | Fig. 3 |
| `figures/fig8_clearance.png` | Fig. 4 |
| `figures/fig5_success.png` | Fig. 5 |
| `figures/fig4_runtime.png` | Fig. 6 |
| `figures/fig6_breakdown.png` | Fig. 7 |
| `figures/fig7_n1.png` | Fig. 8 |
| "Table I: S1 and N1" in `summary.md` | Table III |
| "Per-scenario results" in `summary.md` | Table IV |

## License

MIT. See `LICENSE`.
