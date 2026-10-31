"""
Where things are: the repository, this project's data and runs, and how checkpoints store paths.

The package sits in the CrossGen checkout at py/squadgen/squadgen2d, and everything it reads or writes by default is
in that checkout: the corpus OBJs (data/meshes), the CDF data (py/squadgen/data) and the training runs
(py/squadgen/runs). A checkpoint stores its paths relative to the repository (portable()), so that a run trained on
another machine -- a cloud GPU -- finds its Geom-AE, SQ-VAE and token-set cache here once its run directory has been
copied to the same place under py/squadgen/runs (from_repo()).
"""

import os

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # py/squadgen
REPO = os.path.dirname(os.path.dirname(PROJECT))                        # the CrossGen checkout
DATA = os.path.join(PROJECT, "data")
RUNS = os.path.join(PROJECT, "runs")
# The single-material corpus every default refers to: 35 singlemat + 407 MAMBO domains.
CORPUS = [os.path.join(REPO, "data", "meshes", d) for d in ("singlemat", "mambo")]

DEFAULT_CDF = os.path.join(DATA, "cdf")
DEFAULT_GEOM_AE = os.path.join(RUNS, "geom_ae_2d")
DEFAULT_SQ_VAE = os.path.join(RUNS, "sq_vae_2d")
DEFAULT_SQ_DIFFUSE = os.path.join(RUNS, "sq_diffuse_2d")

# Config entries that are paths. A checkpoint stores them relative to the repository when they are inside it.
PATH_KEYS = ("out", "data", "geom_ae", "vae", "latent_stats", "init", "pack")


def to_repo(path):
    a = os.path.abspath(path)
    return os.path.relpath(a, REPO) if a.startswith(REPO + os.sep) else a


def from_repo(path):
    """A stored path on this machine: a relative one is the repository's; an absolute one that is not here but runs
    through py/squadgen/ (another machine's checkout) is re-rooted at this repository's py/squadgen/."""
    if not os.path.isabs(path):
        return os.path.join(REPO, path)
    marker = os.sep + os.path.relpath(PROJECT, REPO) + os.sep
    if not os.path.exists(path) and marker in path:
        local = os.path.join(PROJECT, path.split(marker, 1)[1])
        if os.path.exists(local):
            return local
    return path


def portable(cfg):
    """cfg with its PATH_KEYS made repository-relative, as checkpoints store it."""
    out = dict(cfg)
    for k in PATH_KEYS:
        v = out.get(k)
        if isinstance(v, str):
            out[k] = to_repo(v)
        elif isinstance(v, (list, tuple)) and v and all(isinstance(x, str) for x in v):
            out[k] = [to_repo(x) for x in v]
    return out
