"""Where ms2exp finds its inputs.

The data directory (--data, else the environment variable MS2EXP_DATA) holds each input under its
logical name, as ms2exp.fetch writes it. The lab root (MS2EXP_LAB_ROOT) is the Langmead lab's
working tree on Rockfish; only the preparation steps that start from lab-internal pools read it,
and they say so when it is unset.
"""
import os
from pathlib import Path


def data_dir(arg=None):
    """The data directory: ARG, else MS2EXP_DATA."""
    value = arg or os.environ.get('MS2EXP_DATA')
    if not value:
        raise SystemExit('Give --data DIR or set MS2EXP_DATA to the directory ms2exp.fetch fills')
    return Path(value).resolve()


def lab_root():
    """The lab's working tree (MS2EXP_LAB_ROOT), for steps that read lab-internal inputs."""
    value = os.environ.get('MS2EXP_LAB_ROOT')
    if not value:
        raise SystemExit('This step reads lab-internal inputs: set MS2EXP_LAB_ROOT (in the Langmead lab, '
                         '/scratch16/blangme2/ms-bakery), or start from the deposited read sets')
    return Path(value)
