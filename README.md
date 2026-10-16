# rottnest-quickstart

Get fault-tolerant resource estimates for **any Cirq circuit** with [Rottnest](https://github.com/QSI-BAQS): Rz rotations, T gates, runtime and space-time volume. It works from any project folder.

## What's here

| File | What it is |
|---|---|
| `rottnest_estimate_template.ipynb` | Notebook to copy into your project. Put in your circuit, click **Run All**. |
| `rottnest_helpers.py` | The bridge between your circuit and Rottnest. The notebook needs it. |
| `Rottnest_Quickstart_Guide.pdf` | 5-page illustrated guide: setup, how it works, troubleshooting, cheat sheet. |

## Quick start

1. Install Rottnest (see [QSI-BAQS](https://github.com/QSI-BAQS)) into a Python 3.11 environment.
2. Make `rottnest_helpers.py` importable: put it next to your notebook, or add its folder to that environment (a `.pth` file in `site-packages`).
3. Copy `rottnest_estimate_template.ipynb` into your project, choose the Rottnest environment as the kernel, and replace the example circuit with yours.
4. Run All.

Or, in plain Python:

```python
from rottnest_helpers import count_cirq_circuit

result = count_cirq_circuit(circuit)   # any cirq.Circuit
print(result["rz_count"], result["t_count"], result["hardware"]["tocks"]["total"])
```

## What the notebook gives you

- Headline numbers: Rz rotations, T gates, T infidelity, runtime (tocks), magic states.
- A table and a six-panel chart of **every** parameter Rottnest reports.
- How each parameter grows with qubit count, if your circuit has a size parameter.
- Saved results: `rottnest_result.json`, `rottnest_parameters.png`, `rottnest_scaling.csv`, `rottnest_scaling.png`.

## Why the helper exists

- **Rottnest runs circuits in separate worker processes,** which can't see notebook variables. The helper passes your circuit as a file and runs each job in a fresh process.
- **Rottnest silently miscounts some gates.** It only understands H, X/Y/Z rotations, Rx/Ry/Rz, T, CNOT and CZ, and treats e.g. `CZ**0.25` as a full CZ. `prepare_circuit` rewrites only the unsupported gates, and the notebook checks the result is the same operation.
- **Measurements make Rottnest hang.** They are removed first; they don't change the counts.

## Notes

- **Architecture:** uses the *Four Stage Superconducting* architecture (t_scheduler). *Rz Counter* currently hangs because of a Rottnest bug.
- **Chip layout:** hardware numbers depend on it. The default is t_scheduler's 20×10 example chip; pass `layout=...` to use another.
- **Tested setup:** macOS with small circuits (up to a few thousand T gates). Large circuits are untested.
- **PDF paths:** the guide describes a setup with Rottnest in `~/Desktop/rottnest`; adjust the paths for yours.
