'''
Helpers for running Rottnest from Python or Jupyter.

    from rottnest_helpers import count_cirq_circuit, run_rottnest

    result = count_cirq_circuit(my_cirq_circuit)          # your own circuit
    result = run_rottnest("Heisenberg", {"height": 2, "width": 2})   # a built-in example

Each run happens in a fresh Python process. Rottnest starts its own pool of
worker processes and only allows one pool per process, so this keeps the
notebook kernel clean and lets you run as many jobs as you like.
'''
import json
import os
import signal
import subprocess
import sys
import tempfile
import time

import cirq

from rottnest.executables.t_rz_executable import T_RZ_RottnestExecutable

DEFAULT_ARCHITECTURE = "Four Stage Superconducting"   # t_scheduler
CIRQ_EXECUTABLE_NAME = "Cirq circuit (JSON file)"
_RESULT_MARKER = "ROTTNEST_RESULT "
# Some Macs have a beta SDK the linker can't read; pin a working one if present.
_SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk"


class CirqJsonCircuit(T_RZ_RottnestExecutable):
    '''
        Rottnest executable that loads any cirq circuit saved with cirq.to_json.
        Rottnest's worker processes re-create executables from their module, so
        the circuit travels as a file path instead of a Python object.
    '''

    @staticmethod
    def get_name():
        return CIRQ_EXECUTABLE_NAME

    @staticmethod
    def get_parameters():
        # Declared directly (like the built-in examples): Rottnest's automatic
        # parameter collection fails on base classes that define none.
        return T_RZ_RottnestExecutable.base_params | {"circuit_file": (str, "")}

    def _generate_circuit(self):
        if getattr(self, "_circ", None) is None:
            self._circ = cirq.read_json(self.circuit_file)
        return self._circ


# Rottnest looks for this module-level list when loading executables.
rottnest_executables = [CirqJsonCircuit]


def count_cirq_circuit(circuit: cirq.Circuit, *, architecture=DEFAULT_ARCHITECTURE,
                       layout=None, timeout=600, verbose=False) -> dict:
    '''
        Run your own cirq circuit through Rottnest and return the counts.
        Measurements are dropped and composite gates (e.g. cirq.qft) are
        expanded first: Rottnest hangs on both, and neither affects the counts.
    '''
    circuit = prepare_circuit(circuit)
    fd, path = tempfile.mkstemp(suffix=".json", prefix="rottnest_circuit_")
    os.close(fd)
    try:
        cirq.to_json(circuit, path)
        return run_rottnest(
            CIRQ_EXECUTABLE_NAME,
            {"circuit_file": path},
            architecture=architecture,
            layout=layout,
            extra_modules=["rottnest_helpers"],
            timeout=timeout,
            verbose=verbose,
        )
    finally:
        os.remove(path)


def _rottnest_supports(op: cirq.Operation) -> bool:
    '''
        Gates Rottnest's cirq parser handles correctly (rottnest/monkey_patchers/
        cirq_patcher.py). It treats CZ**t / CNOT**t as full CZ / CNOT whatever t
        is, so only whole powers count as supported.
    '''
    gate = op.gate
    if isinstance(gate, (cirq.CZPowGate, cirq.CXPowGate)):
        return gate.global_shift == 0 and float(gate.exponent) % 2 == 1
    if isinstance(gate, cirq.HPowGate):
        return float(gate.exponent) % 2 == 1
    return isinstance(gate, (cirq.XPowGate, cirq.YPowGate, cirq.ZPowGate))


def _compact_form(op: cirq.Operation) -> list | None:
    '''
        Textbook decompositions for common gates, which use far fewer rotations
        than cirq's general-purpose synthesis.
    '''
    gate = op.gate
    if isinstance(gate, cirq.CZPowGate) and gate.global_shift == 0:
        # Controlled phase: Z^(t/2) on both, then CNOT . Z^(-t/2) . CNOT
        c, tgt = op.qubits
        h = float(gate.exponent) / 2
        return [(cirq.Z ** h).on(c), (cirq.Z ** h).on(tgt), cirq.CNOT(c, tgt),
                (cirq.Z ** -h).on(tgt), cirq.CNOT(c, tgt)]
    if isinstance(gate, cirq.CXPowGate) and gate.global_shift == 0:
        c, tgt = op.qubits
        return [cirq.H(tgt), *_compact_form((cirq.CZ ** gate.exponent).on(c, tgt)),
                cirq.H(tgt)]
    if isinstance(gate, cirq.SwapPowGate) and float(gate.exponent) % 2 == 1:
        a, b = op.qubits
        return [cirq.CNOT(a, b), cirq.CNOT(b, a), cirq.CNOT(a, b)]
    return None


def prepare_circuit(circuit: cirq.Circuit) -> cirq.Circuit:
    '''
        Make a circuit safe for Rottnest without changing what it does:
        drop measurements, and rewrite every gate Rottnest doesn't support
        (controlled-phase rotations, SWAP, iSWAP, FSim, composite gates such as
        cirq.qft, ...) into full CZ gates plus single-qubit X/Y/Z rotations.
        Supported gates are passed through untouched.
    '''
    ops = []
    keep = lambda o: cirq.is_measurement(o) or _rottnest_supports(o)
    for op in cirq.decompose(circuit, keep=keep, on_stuck_raise=None):
        if cirq.is_measurement(op) or not op.qubits:
            continue
        if _rottnest_supports(op):
            ops.append(op)
            continue
        compact = _compact_form(op)
        if compact is not None:
            ops.extend(compact)
            continue
        native = cirq.optimize_for_target_gateset(
            cirq.Circuit(op), gateset=cirq.CZTargetGateset()
        )
        for sub in native.all_operations():
            if not sub.qubits:
                continue
            if _rottnest_supports(sub):
                ops.append(sub)
            elif len(sub.qubits) == 1:
                q = sub.qubits[0]
                for pauli, half_turns in cirq.single_qubit_matrix_to_pauli_rotations(
                        cirq.unitary(sub), atol=1e-9):
                    ops.append((pauli ** half_turns).on(q))
            else:
                raise ValueError(f"Can't convert {sub} for Rottnest")
    return cirq.Circuit(ops)


def run_rottnest(executable: str, params: dict | None = None, *,
                 architecture=DEFAULT_ARCHITECTURE, layout=None,
                 extra_modules=(), timeout=600, verbose=False) -> dict:
    '''
        Run one Rottnest job (same pipeline as the web app's Run button).

        executable    : name shown by list_options(), e.g. "Heisenberg"
        params        : executable parameters, e.g. {"height": 2, "width": 2}
        architecture  : "Four Stage Superconducting" (t_scheduler).
                        "Rz Counter" currently hangs because of a Rottnest bug.
        layout        : chip layout dict; defaults to t_scheduler's 20x10 example
        Returns a dict with rz_count, t_count, t_infidelity, hardware results
        (space-time volumes, tocks, T sources) and the run time.
    '''
    job = {
        "executable": executable,
        "params": params or {},
        "architecture": architecture,
        "layout": layout,
        "extra_modules": list(extra_modules),
    }
    env = dict(os.environ)
    if os.path.isdir(_SDKROOT):   # macOS workaround; ignored elsewhere
        env["SDKROOT"] = _SDKROOT
    here = os.path.dirname(os.path.abspath(__file__))
    env["PYTHONPATH"] = here + os.pathsep + env.get("PYTHONPATH", "")

    # Rottnest leaves gridsynth helper processes behind, so the job runs in its
    # own process group (cleaned up below) and writes to a file, not a pipe.
    with tempfile.TemporaryFile(mode="w+") as log:
        proc = subprocess.Popen(
            [sys.executable, os.path.join(here, "rottnest_helpers.py"), json.dumps(job)],
            stdout=log, stderr=subprocess.STDOUT, env=env, cwd=here,
            start_new_session=True,
        )
        try:
            proc.wait(timeout=timeout)
        finally:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        log.seek(0)
        output = log.read()

    if verbose:
        print(output)
    for line in output.splitlines():
        if line.startswith(_RESULT_MARKER):
            return json.loads(line[len(_RESULT_MARKER):])
    raise RuntimeError(
        f"Rottnest run failed (exit code {proc.returncode}).\n"
        + "\n".join(output.splitlines()[-30:])
    )


def list_options() -> dict:
    '''Executables and architectures Rottnest can currently see.'''
    from rottnest.plugins import executables, architectures
    executables.load_modules_from_strings("rottnest_helpers")
    return {
        "executables": executables.get_executable_names(),
        "architectures": architectures.get_architecture_names(),
    }


def _run_job_in_this_process(job: dict) -> dict:
    from rottnest.plugins import executables, architectures
    from rottnest.procedures import preprocess_and_execute
    from rottnest.process_pool import get_pool
    from rottnest.compute_units.layout_proxy import LayoutProxy
    from t_scheduler.region_builder.json_to_region import example

    if job["extra_modules"]:
        executables.load_modules_from_strings(*job["extra_modules"])
    if not executables.set_current_executable(job["executable"]):
        raise ValueError(f"Unknown executable {job['executable']!r}")
    current = executables.get_executable_params()
    executables.set_executable_params(**(current | job["params"]))
    architectures.set_current_architecture(job["architecture"])
    LayoutProxy.add_layout_with_id(0, job["layout"] or example)

    t0 = time.time()
    procedure = preprocess_and_execute.PreprocAndExecuteProcedure(reporting=False)
    try:
        procedure.execute()
        while not procedure.complete():
            procedure.poll()

        pre = procedure.preprocessor
        hardware = procedure.pool_procedure.get_results.get_results()
        if hasattr(hardware, "to_runchart"):
            hardware = hardware.to_runchart()
        return {
            "executable": job["executable"],
            "architecture": job["architecture"],
            "rz_count": pre.get_rz_count(),
            "t_count": pre.get_t_count(),
            "t_infidelity": pre.get_t_infidelity(),
            "hardware": hardware,
            "seconds": round(time.time() - t0, 2),
        }
    finally:
        # Without this the worker processes keep the job alive forever.
        get_pool().terminate()


if __name__ == "__main__":
    try:
        _result = _run_job_in_this_process(json.loads(sys.argv[1]))
        print(_RESULT_MARKER + json.dumps(_result, default=str), flush=True)
        _code = 0
    except Exception:
        import traceback
        traceback.print_exc()
        _code = 1
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(_code)   # don't wait on Rottnest's leftover helper threads
