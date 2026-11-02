"""Command execution under /usr/bin/time with the benchmark environment, and atomic JSON writes."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

# The repository root goes on PYTHONPATH so that commands run as python3 -m ms2exp... from any directory.
REPO = Path(__file__).resolve().parents[2]
ENV = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
       'PYTHONPATH': os.pathsep.join(p for p in (str(REPO), os.environ.get('PYTHONPATH')) if p)}


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def execute(argv, folder, label, timeout=3600, env=None):
    """Block in wait(), not Python's timeout polling loop, for precise wall time."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    argv = list(map(str, argv))
    out, err, metrics = [folder / (label + ext) for ext in ['.stdout', '.stderr', '.time']]
    command = (['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(metrics)] + argv
               if sys.platform.startswith('linux') else argv)
    expired = threading.Event()
    with out.open('w') as stdout, err.open('w') as stderr:
        started = time.perf_counter()
        process = subprocess.Popen(command, cwd=folder, env=env or ENV,
                                   stdout=stdout, stderr=stderr, start_new_session=True)
        def kill():
            expired.set()
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        timer = threading.Timer(timeout, kill)
        timer.daemon = True
        timer.start()
        try:
            rc = process.wait()
            elapsed = time.perf_counter() - started
        finally:
            timer.cancel()
    row = dict(argv=argv, cwd=str(folder), returncode=rc, timed_out=expired.is_set(),
               wall_seconds=elapsed, stdout=str(out), stderr=str(err))
    if metrics.exists():
        lines = metrics.read_text().splitlines()
        fields = lines[-1].split() if lines else []
        if len(fields) == 4:
            row.update(native_wall_seconds=float(fields[0]), user_seconds=float(fields[1]),
                       system_seconds=float(fields[2]), peak_rss_kib=int(fields[3]))
    save(folder / (label + '.command.json'), row)
    if rc or expired.is_set():
        raise RuntimeError(f'{label}: exit {rc}, timeout={expired.is_set()}; see {err}')
    return row


def atomic_save(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)
