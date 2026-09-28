#!/usr/bin/env python3
"""Predict, for each block-decomposition method, whether it will produce a valid decomposition of a
mesh and how good that decomposition will be, from the mesh's boundary features and the result of
one cheap method run first (the probe); rank the methods by the expected result.

DATA  dataset.csv from build_dataset.py (CSV, or Parquet if the file name ends in .parquet):
        mesh_id, group, source, area, [lam_1..lam_K], feat_*, probe_<PROBE>_*, <method>_<column>
  inputs    the column groups in INPUTS that the data has: feat_*, Mesh.boundary_features()
            (corners, holes, T3's bound; area-free); probe_<PROBE>_*, one run of the PROBE method
            (below); and, if INPUTS asks for it, lam_1..lam_K, the spectrum in build_dataset.py's
            "weyl_ratio" normalisation (area-free, divided by Weyl's line). Each column is
            transformed as its range calls for -- log if positive, log1p if a count, asinh if it
            can go negative -- and standardised. A row with a missing input is skipped, except for
            the probe's metrics and block count, which are missing when it was not valid or raised:
            they count as 0, and its _valid column says which.
  outputs   per method, <method>_valid, the fraction of its runs that produced a decomposition, and
            the METRICS of a valid one (each in [0, 1], higher better: regularity, angle_quality,
            chord_quality). The two are learnt by separate heads: a validity logit (binary cross-
            entropy) and the quality given validity (MSE on runs that were valid). Data from before
            build_dataset.py wrote _valid, where a failed method scored 0 on every metric, is read
            the same way: all-zero is not valid. An empty value is masked out.
  utility   of one method on one mesh: p * U(quality) + (1 - p) * U(failure), p the probability it
            is valid, U the METRICS weighted as below and a failure counting as their failed values;
            with FALLBACK set, a failed run is taken to be detected (BlockDecomposition.valid) and
            replaced by that method's run, so U(failure) is the fallback's own expected utility.
  probe     PROBE runs before the selector is asked, and its result is an input. The model's own
            prediction for it is replaced by that result, and the model is meant to be used as:
            run PROBE; if it was valid and ranks first, keep it; otherwise run the best-ranked other
            method and keep the better of the two (FALLBACK's run if neither is valid). The report
            scores exactly that ("<probe> first"), with its seconds per mesh.
  group     remeshes / faces of one part share a group (build_dataset.py sets it), so a part never
            lands on both sides of a train/test split. The report also cross-validates by shape:
            meshes whose boundary features are identical form one group, which is the test for a
            shape the model has not seen -- identical faces recur across parts, so the part groups
            alone let a model be tested on shapes it trained on. The shape groups are geometry only
            (the spectrum where there are no features, never the probe): two copies of one face are
            one shape even where the probe ran differently on them.

USAGE
  python train_classifier.py --data dataset.csv [--onnx selector.onnx] [--plot live|save|off]
                             [--cv both|part|shape]
      cross-validated report, then trains on all meshes and saves selector.pt. While it trains, a
      window shows each net's training and validation loss per epoch (--plot live, the default), and
      each stage's curves are saved to --plot-dir as loss_<stage>.png (save: the PNGs only).
  python train_classifier.py --predict --data new.csv
      ranks the methods for new meshes (needs mesh_id and the model's input columns, the probe's
      result among them) -> predictions.csv, with the probe rule's action for each mesh

All preprocessing is inside the saved model: it takes the raw input columns and returns the
predicted quality of a valid layout (method, metric), the probability of one per method, and a
utility per method. Needs numpy, pandas and torch; Parquet needs pyarrow, --onnx needs onnx and
onnxscript, --plot needs matplotlib.
"""
import argparse
import copy
import logging
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

# ------------------------------------------------------------------------------------------------
# EDIT to match your data.   metric name: (better, log_transform, utility_weight, failed_value)
#   better          "higher" or "lower"
#   log_transform   True for positive, skewed metrics (errors, runtimes)
#   utility_weight  importance in the ranking. Metrics are standardized over valid runs, so weights
#                   are in standard deviations: 1.0 vs 0.5 means one std of the first is worth two
#                   std of the second. Changing weights later needs no retraining (--predict uses
#                   the current values).
#   failed_value    what a run that produced no valid decomposition counts as on this metric, in
#                   its own units (for the utility only; it is never trained on). 0 is the worst
#                   value of the three [0, 1] metrics: no layout is the worst layout.
# angle_quality has no weight: it is read on the method's own polylines before smoothing moves them
# -- the metrics note's B4, "report it, do not rank on it" -- and it is the least repeatable of the
# three across copies of one shape. It is still predicted, as a second task for the quality head.
# The methods are every <method>_<metric> prefix in the data, in column order.
METRICS = {
    "regularity": ("higher", False, 1.0, 0.0),
    "angle_quality": ("higher", False, 0.0, 0.0),
    "chord_quality": ("higher", False, 1.0, 0.0),
}
# Input column groups, used when the data has them: "feat" is every feat_*, "probe" every
# probe_<PROBE>_*, "lam" is lam_1..lam_K. Measured on the 2026-09-28 data (442 meshes; 5-fold CV by
# shape, mean over 3-4 fold assignments, each moving it by ~0.01; lower better): regret with the
# fallback, feat 0.186, feat + lam 0.219, lam 0.265, feat + probe 0.110; with PROBE's rule, feat +
# probe 0.086 at 3.6 s a mesh (the rule with no model, zipline + atlas, 0.121 at 6.4 s). The spectrum
# also costs a median 11 s a mesh to compute, more than any one method, so it is not used; add "lam"
# to use it.
INPUTS = ("feat", "probe")
# A method whose run replaces one that failed (None: a failure is final). ATLAS never failed on the
# 2026-09-26 or 2026-09-28 data, and ranking with it cut the regret to half of ranking without.
FALLBACK = "atlas"
# The method that runs before the selector is asked, whose result is then an input (None: none).
# ZIPLINE takes ~0.1 s, and as a cross-field method its block count and quality say whether MERIDIAN
# and TORSION will split a face into many blocks, which the boundary features cannot. Its prediction
# is replaced by its result, and the model is for the rule: run PROBE; keep it if it was valid and
# ranks first; else run the best-ranked other method and keep the better of the two.
PROBE = "zipline"
HIDDEN = [64, 64]    # MLP hidden layer widths ([] = linear model)
DROPOUT = 0.1
LR, WEIGHT_DECAY, BATCH = 1e-3, 1e-2, 32
MAX_EPOCHS, PATIENCE = 2000, 150   # early stopping on a held-out ~15% of groups
PLOT_EVERY = 10      # epochs between redraws of the live loss plot
REL_WEIGHT = 1.0     # extra loss weight on between-method differences (they decide the ranking)
N_SEEDS, N_FOLDS = 5, 5            # ensemble size, cross-validation folds
SAME_SHAPE = 1e-3    # RMS difference of standardized inputs below which two meshes are one shape
# ------------------------------------------------------------------------------------------------

torch.set_num_threads(1)           # tiny networks train faster on one thread

LOG, LOG1P, ASINH = 1, 2, 3


# ---------------------------------------------------------------------------------------- model
class Features(nn.Module):
    """Raw input row -> transformed, column by column: log for a positive column (the weyl_ratio
    spectrum sits near 1 for large k but lam_1 runs from ~1.4 to ~18, so a log makes a relative
    difference the same size at every k), log1p for a non-negative one (corner counts), asinh for
    one that can go negative (the Euler characteristic). The kinds come from the training data and
    are saved with the model. The row must not be sorted: lambda_k / k is not monotone in k."""

    def __init__(self, kinds):
        super().__init__()
        self.register_buffer("kind", torch.tensor(kinds, dtype=torch.int64))

    def forward(self, x):
        k = self.kind
        # asinh written out: the TorchScript ONNX exporter has no symbol for aten::asinh.
        asinh = torch.sign(x) * torch.log(x.abs() + torch.sqrt(x * x + 1.0))
        return torch.where(k == LOG, torch.log(x.clamp_min(1e-30)),
                           torch.where(k == LOG1P, torch.log1p(x.clamp_min(0.0)), asinh))


class Net(nn.Module):
    """Features -> standardize -> MLP -> per method a validity logit, and standardized quality
    scores (batch, method, metric), higher = better."""

    def __init__(self, cfg):
        super().__init__()
        self.feats = Features(cfg["input_kinds"])
        d = len(cfg["input_cols"])
        self.register_buffer("mu", torch.zeros(d))
        self.register_buffer("sd", torch.ones(d))
        layers = []
        for h in cfg["hidden"]:
            layers += [nn.Linear(d, h), nn.GELU(), nn.Dropout(cfg["dropout"])]
            d = h
        self.A, self.M = len(cfg["algorithms"]), len(cfg["metrics"])
        self.mlp = nn.Sequential(*layers, nn.Linear(d, self.A * (1 + self.M)))

    def forward(self, x):
        out = self.mlp((self.feats(x) - self.mu) / self.sd)
        return out[:, :self.A], out[:, self.A:].reshape(-1, self.A, self.M)


class Selector(nn.Module):
    """Ensemble of Nets + target transform. Input: raw input columns, float32 (batch, D); the
    probe's metrics and block count may be NaN. Output: metrics (batch, method, metric), the
    predicted quality of a valid decomposition in original units -- flattened, the <method>_<metric>
    columns in their order; valid (batch, method), the probability of a valid one; utility (batch,
    method), the expected utility. Rank methods by descending utility. The probe's own entries are
    its result, and its utility is -inf when it was not valid: it is no longer an option."""

    def __init__(self, cfg):
        super().__init__()
        self.nets = nn.ModuleList(Net(cfg) for _ in range(cfg["n_members"]))
        a, m = len(cfg["algorithms"]), len(cfg["metrics"])
        self.register_buffer("y_mu", torch.zeros(m))
        self.register_buffer("y_sd", torch.ones(m))
        self.register_buffer("sign", torch.tensor(cfg["sign"], dtype=torch.float32))
        self.register_buffer("w", torch.tensor(cfg["weights"], dtype=torch.float32))
        self.register_buffer("is_log", torch.tensor(cfg["is_log"]))
        self.register_buffer("t_fail", torch.tensor(cfg["t_fail"], dtype=torch.float32))
        self.register_buffer("fallback", torch.tensor(cfg["fallback"], dtype=torch.int64))
        # Rebuilt from cfg rather than saved, so a model from before the probe still loads (no
        # probe, nothing filled). probe_index: the input columns of its _valid, then its metrics.
        probe = cfg.get("probe", -1)
        self.register_buffer("fill", torch.tensor(cfg.get("fill", [False] * len(cfg["input_cols"]))),
                             persistent=False)
        self.register_buffer("is_probe", torch.arange(a) == probe, persistent=False)
        self.register_buffer("probe_index", torch.tensor(cfg.get("probe_index", [0] * (1 + m)), dtype=torch.int64),
                             persistent=False)

    def transformed(self, x):
        """Probability of a valid run p, standardized scores z, and metrics in transformed (log
        where configured) units t, the probe's taken from its result in x."""
        x = torch.where(self.fill & torch.isnan(x), torch.zeros_like(x), x)
        outs = [net(x) for net in self.nets]
        p = torch.stack([torch.sigmoid(o[0]) for o in outs]).mean(0)
        z = torch.stack([o[1] for o in outs]).mean(0)
        # Tensor ops, not an `if`, so that an exported graph keeps the substitution.
        pv = torch.index_select(x, 1, self.probe_index[:1])
        pt = torch.index_select(x, 1, self.probe_index[1:])
        pz = self.sign * (torch.where(self.is_log, torch.log(pt.clamp_min(1e-30)), pt) - self.y_mu) / self.y_sd
        p = torch.where(self.is_probe, pv, p)
        z = torch.where(self.is_probe[:, None], pz[:, None, :], z)
        return p, z, self.sign * z * self.y_sd + self.y_mu

    def forward(self, x):
        p, z, t = self.transformed(x)
        u_valid = (z * self.w).sum(-1)
        u_fail = (self.w * self.sign * (self.t_fail - self.y_mu) / self.y_sd).sum()
        u = p * u_valid + (1 - p) * u_fail
        # With a fallback, a failed run is worth the fallback method's own expected utility.
        # Selected with tensor ops rather than an `if`, so that an exported graph keeps the switch.
        u_fb = torch.index_select(u, 1, self.fallback.clamp_min(0).reshape(1))
        u = torch.where(self.fallback >= 0, p * u_valid + (1 - p) * u_fb, u)
        u = torch.where(self.is_probe & (p < 0.5), torch.full_like(u, float("-inf")), u)
        return torch.where(self.is_log, torch.exp(t), t), p, u


# ----------------------------------------------------------------------------------------- data
def read_table(path):
    return pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path)


def input_columns(D, groups):
    """The columns of each input group the table has, in order: lam_k by k, feat_* and PROBE's
    probe_* as they come."""
    cols = []
    if "lam" in groups:
        cols += sorted((c for c in D.columns if re.fullmatch(r"lam_\d+", str(c))), key=lambda c: int(c[4:]))
    if "feat" in groups:
        cols += [c for c in D.columns if str(c).startswith("feat_")]
    if "probe" in groups and PROBE:
        cols += [c for c in D.columns if str(c).startswith(f"probe_{PROBE}_")]
    return cols


def fillable(cols):
    """Per column, whether it may be missing: a probe's metrics and block count, when it was not
    valid or raised. They count as 0 -- the metrics' failed value, and no blocks -- and the probe's
    _valid column, which may not be missing, says why."""
    return [str(c).startswith("probe_") and not str(c).endswith("_valid") for c in cols]


def load_data(path, input_cols=None):
    """The table and its input columns (INPUTS' groups, or the model's own when predicting). Rows
    with a missing input (other than a fillable() one) are dropped: there is nothing to predict
    from."""
    D = read_table(path)
    if "mesh_id" not in D.columns:
        sys.exit(f"{path}: needs a mesh_id column")
    D["mesh_id"] = D["mesh_id"].astype(str)
    if D["mesh_id"].duplicated().any():
        sys.exit(f"{path}: duplicate mesh_id values")
    if input_cols is None:
        input_cols = input_columns(D, INPUTS)
        if not input_cols:
            sys.exit(f"{path}: no input columns of the groups {INPUTS} (lam_1..lam_K, feat_*, probe_{PROBE}_*)")
        missing = [g for g in INPUTS if not input_columns(D, (g,))]
        if missing:
            print(f"note: the data has no {' or '.join(missing)} columns; training on the rest"
                  + (f" (build_dataset.py writes probe_{PROBE}_*; rerun it on its cache, which only rewrites "
                     "the CSVs)" if "probe" in missing and PROBE else ""))
    else:
        missing = [c for c in input_cols if c not in D.columns]
        if missing:
            sys.exit(f"{path}: missing columns the model was trained on: {missing[:5]}")
    need = [c for c, f in zip(input_cols, fillable(input_cols)) if not f]
    bad = D[need].isna().any(axis=1).to_numpy()
    if bad.any():
        print(f"note: {bad.sum()} meshes have a missing input and are skipped")
        D = D[~bad].reset_index(drop=True)
    return D, input_cols


def input_kinds(X):
    """Each column's transform, from its range in the training data."""
    return [LOG if (x > 0).all() else LOG1P if (x >= 0).all() else ASINH for x in X.T]


def load_targets(D, metrics):
    """Methods (the <method>_<metric> prefixes, in column order, not counting input columns such as
    probe_zipline_regularity); Y (mesh, method, metric), the quality of a valid run, NaN where
    there is none; V (mesh, method), the fraction of runs that were valid, NaN where unknown."""
    algs = []
    for c in map(str, D.columns):
        if c.startswith(("lam_", "feat_", "probe_")):
            continue
        for m in metrics:
            if c.endswith("_" + m) and c[:-len(m) - 1] not in algs:
                algs.append(c[:-len(m) - 1])
    cols = [f"{a}_{m}" for a in algs for m in metrics]
    missing = [c for c in cols if c not in D.columns]
    if not algs or missing:
        sys.exit(f"need a <method>_<metric> column for every method and each of {metrics}"
                 + (f"; missing {missing[:5]}" if missing else "") + "; edit METRICS")
    Y = D[cols].to_numpy(float).reshape(len(D), len(algs), len(metrics))
    if all(f"{a}_valid" in D.columns for a in algs):
        V = D[[f"{a}_valid" for a in algs]].to_numpy(float)
    else:   # before build_dataset.py wrote _valid: a failed method scored 0 on every metric
        V = np.where(np.isnan(Y).all(2), np.nan, (~(Y == 0).all(2)).astype(float))
        print("note: no <method>_valid columns; a method that scored 0 on every metric is taken as not valid")
    Y = np.where((V > 0)[..., None], Y, np.nan)
    return algs, Y, V


def load_seconds(D, algs):
    """Seconds per run (mesh, method) from the <method>_seconds columns, or None without them. For
    the report's cost of each way of choosing; never trained on."""
    cols = [f"{a}_seconds" for a in algs]
    return D[cols].to_numpy(float) if all(c in D.columns for c in cols) else None


def transform(Y, cfg):
    T = Y.copy()
    for j, name in enumerate(cfg["metrics"]):
        if cfg["is_log"][j]:
            if (T[..., j] <= 0).any():
                sys.exit(f"metric '{name}' has values <= 0, so it can't be log-transformed")
            T[..., j] = np.log(T[..., j])
    return T


def failed_values(cfg):
    """Each metric's failed_value in transformed units."""
    out = []
    for name, is_log in zip(cfg["metrics"], cfg["is_log"]):
        v = METRICS[name][3] if name in METRICS else 0.0
        if is_log and v <= 0:
            sys.exit(f"metric '{name}' is log-transformed, so its failed_value must be positive")
        out.append(float(np.log(v)) if is_log else float(v))
    return out


def same_shape_groups(Z):
    """Connected components of 'the same inputs': rows of the standardized input matrix Z closer
    than SAME_SHAPE (RMS over columns). In blocks, so the distance matrix is never held whole."""
    n, d = Z.shape
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    sq = (Z * Z).sum(1)
    for lo in range(0, n, 1024):
        hi = min(n, lo + 1024)
        d2 = (sq[lo:hi, None] + sq[None, :] - 2.0 * Z[lo:hi] @ Z.T) / d
        for i, j in zip(*np.nonzero(d2 < SAME_SHAPE ** 2)):
            a, b = find(lo + i), find(j)
            if a != b:
                parent[a] = b
    return np.array([f"shape{find(i)}" for i in range(n)])


def shape_groups(D):
    """Meshes of one shape, by same_shape_groups on the table's boundary features (its spectrum
    where it has none; each mesh on its own where it has neither). Geometry only, whatever INPUTS
    is: a probe's result is not the shape, and two copies of one face whose probe runs differ are
    still one shape, which grouping by every input would split between training and test."""
    cols = input_columns(D, ("feat",)) or input_columns(D, ("lam",))
    groups = D["mesh_id"].to_numpy(object).copy()
    G = D[cols].to_numpy(np.float64)
    ok = ~np.isnan(G).any(1) if cols else np.zeros(len(D), bool)
    if ok.sum() > 1:
        with torch.no_grad():
            Fx = Features(input_kinds(G[ok]))(torch.tensor(G[ok])).numpy()
        groups[ok] = same_shape_groups((Fx - Fx.mean(0)) / (Fx.std(0) + 1e-6))
    return groups.astype(str)


# ------------------------------------------------------------------------------------- plotting
class LossPlot:
    """Training and validation loss per epoch, drawn while the nets train. One stage (a cross-
    validated model, or the final fit) is one figure: the net in training is in colour, the ones
    before it in the stage fade to grey, and the stage is saved as loss_<stage>.png when it ends.
    Both losses are loss_fn in eval mode (no dropout), so the two curves are measured alike and the
    gap between them is overfitting, not dropout noise."""

    def __init__(self, out_dir, live):
        import matplotlib
        if not live:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        self.plt = plt
        self.live = live and matplotlib.get_backend().lower() != "agg"
        if live and not self.live:
            print("note: no interactive matplotlib backend; loss plots are only saved")
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.fig, self.ax = plt.subplots(figsize=(9, 5.5))
        if self.live:
            plt.ion()
            plt.show(block=False)

    def begin_stage(self, name):
        self.stage = name
        self.ax.clear()
        self.ax.set_yscale("log")
        self.ax.set_xlabel("epoch")
        self.ax.set_ylabel("loss (validity BCE + standardised quality MSE + between-method term)")
        self.ax.grid(True, which="both", alpha=0.3)
        self.cur = None

    def begin_net(self, label):
        if self.cur:     # the previous net stays on the plot, greyed out and out of the legend
            for line in self.cur:
                line.set(color="0.75", alpha=0.6, label="_")
        self.label, self.tr, self.va = label, [], []
        self.cur = [self.ax.plot([], [], color="C0", lw=1.5, label="train")[0],
                    self.ax.plot([], [], color="C1", lw=1.5, label="validation")[0],
                    self.ax.axvline(0, color="C2", ls="--", lw=1, label="best (kept) epoch")]
        self.ax.legend(loc="upper right")

    def epoch(self, train_loss, val_loss, best_epoch):
        self.tr.append(train_loss)
        self.va.append(val_loss)
        self.best_epoch = best_epoch
        if len(self.tr) % PLOT_EVERY == 0:
            self.redraw()

    def redraw(self):
        x = np.arange(1, len(self.tr) + 1)
        tr, va, best = self.cur
        tr.set_data(x, self.tr)
        va.set_data(x, self.va)
        best.set_xdata([self.best_epoch + 1] * 2)
        self.ax.relim()
        self.ax.autoscale_view()
        self.ax.set_title(f"{self.stage} | {self.label} | epoch {len(self.tr)}, best {self.best_epoch + 1} "
                          f"(validation {self.va[self.best_epoch]:.4f})", fontsize=10)
        if self.live:
            # Not plt.pause: that re-shows the window and on macOS pulls it to the front each time.
            self.fig.canvas.draw_idle()
            self.fig.canvas.start_event_loop(0.001)

    def end_net(self):
        self.redraw()

    def end_stage(self):
        path = self.out_dir / ("loss_" + re.sub(r"[^A-Za-z0-9]+", "_", self.stage).strip("_") + ".png")
        self.fig.savefig(path, dpi=130, bbox_inches="tight")
        print(f"  loss curves -> {path}")


# ------------------------------------------------------------------------------------- training
def group_folds(groups, k, seed):
    """Assign whole groups (not single meshes) to k folds."""
    uniq = np.unique(groups)
    np.random.default_rng(seed).shuffle(uniq)
    fold = {g: i % k for i, g in enumerate(uniq)}
    return np.array([fold[g] for g in groups])


def loss_fn(out, v, mv, z, mq):
    """Binary cross-entropy of the validity logits against the fraction of valid runs v (mask mv),
    plus MSE on the standardized quality of valid runs z (mask mq) and REL_WEIGHT * MSE on its
    between-method differences. Kept apart because a failure is an event, not a quality: as a 0 on
    every metric it was most of each metric's variance, and a squared error spent the net on it."""
    logit, q = out
    bce = (F.binary_cross_entropy_with_logits(logit, v, reduction="none") * mv).sum() / mv.sum().clamp_min(1)
    e = (q - z) * mq
    e_rel = (e - e.sum(1, keepdim=True) / mq.sum(1, keepdim=True).clamp_min(1)) * mq
    return bce + (e.pow(2).sum() + REL_WEIGHT * e_rel.pow(2).sum()) / mq.sum().clamp_min(1)


def train_net(net, X, targets, groups, seed, plot=None, label=""):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    with torch.no_grad():
        Fx = net.feats(X)
        net.mu.copy_(Fx.mean(0))
        net.sd.copy_(Fx.std(0).clamp_min(1e-6))
    val = group_folds(groups, 7, seed) == 0
    tr, va = np.flatnonzero(~val), torch.from_numpy(np.flatnonzero(val))
    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    tr_t = torch.from_numpy(tr)
    at = lambda idx: [t[idx] for t in targets]  # noqa: E731
    best, best_state, best_epoch, bad = float("inf"), None, 0, 0
    if plot:
        plot.begin_net(label)
    for epoch in range(MAX_EPOCHS):
        net.train()
        for b in np.array_split(rng.permutation(tr), max(1, len(tr) // BATCH)):
            b = torch.from_numpy(b)
            loss = loss_fn(net(X[b]), *at(b))
            opt.zero_grad()
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            v = loss_fn(net(X[va]), *at(va)).item()
            t = loss_fn(net(X[tr_t]), *at(tr_t)).item() if plot else None
        if v < best - 1e-5:
            best, best_state, best_epoch, bad = v, copy.deepcopy(net.state_dict()), epoch, 0
        else:
            bad += 1
        if plot:
            plot.epoch(t, v, best_epoch)
        if bad >= PATIENCE:
            break
    if plot:
        plot.end_net()
    net.load_state_dict(best_state)


def fit_selector(cfg, X, T, V, groups, seed=0, plot=None, label=""):
    """Train an ensemble on transformed quality targets T (mesh, method, metric; NaN where no valid
    run) and validity fractions V (mesh, method; NaN = unknown)."""
    mu, sd = np.nanmean(T, axis=(0, 1)), np.nanstd(T, axis=(0, 1)) + 1e-12
    Z = np.array(cfg["sign"]) * (T - mu) / sd   # pooled over methods, so they stay comparable
    torch.manual_seed(seed)    # the initial weights: the same call trains the same model
    sel = Selector(cfg)
    sel.y_mu.copy_(torch.tensor(mu))
    sel.y_sd.copy_(torch.tensor(sd))
    Xt = torch.tensor(X, dtype=torch.float32)
    targets = [torch.tensor(np.nan_to_num(V), dtype=torch.float32),
               torch.tensor(~np.isnan(V), dtype=torch.float32),
               torch.tensor(np.nan_to_num(Z), dtype=torch.float32),
               torch.tensor(~np.isnan(Z), dtype=torch.float32)]
    for i, net in enumerate(sel.nets):
        train_net(net, Xt, targets, groups, seed + i, plot, f"{label}member {i + 1}/{len(sel.nets)}")
    return sel.eval()


def cv_predictions(cfg, X, T, V, groups, plot=None, label=""):
    """Out-of-fold P(valid) (mesh, method) and quality in transformed units (mesh, method, metric)."""
    folds = group_folds(groups, N_FOLDS, seed=12345)
    P = np.full_like(V, np.nan)
    Tp = np.full(T.shape, np.nan)
    for f in range(N_FOLDS):
        tr = folds != f
        sel = fit_selector(cfg, X[tr], T[tr], V[tr], groups[tr], seed=100 * f, plot=plot,
                           label=f"{label}fold {f + 1}/{N_FOLDS}, ")
        with torch.no_grad():
            p, _, t = sel.transformed(torch.tensor(X[~tr], dtype=torch.float32))
        P[~tr], Tp[~tr] = p.numpy(), t.numpy()
    return P, Tp


# ------------------------------------------------------------------------------------ evaluation
def utilities(P, T, cfg, mu, sd, fallback):
    """Expected utility (mesh, method) of choosing each method, from the probability of a valid run
    P and the quality of one T (transformed units), with the metrics standardized by mu, sd."""
    k = np.array(cfg["sign"]) * np.array(cfg["weights"]) / sd
    u_valid = np.nan_to_num(((T - mu) * k).sum(-1))
    u_fail = float(((np.array(cfg["t_fail"]) - mu) * k).sum())
    U = P * u_valid + (1 - P) * u_fail
    if fallback is not None:
        U = P * u_valid + (1 - P) * U[:, [fallback]]
    return U


def avg_rank(u):
    """Ranks along the last axis, ties sharing their average rank."""
    r = u.argsort(1).argsort(1).astype(float)
    for i in range(len(u)):
        for val in np.unique(u[i]):
            same = u[i] == val
            if same.sum() > 1:
                r[i, same] = r[i, same].mean()
    return r


def auc(score, label):
    pos, neg = score[label], score[~label]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    r = avg_rank(np.concatenate([pos, neg])[None, :])[0] + 1
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def selection(U, Up, tol=1e-9):
    """Regret and friends of picking argmax Up where U is what each method is actually worth. A
    pick is right when it is within tol of the best: on the many meshes where methods tie, any of
    the tied ones is."""
    i, pick, best = np.arange(len(U)), Up.argmax(1), U.max(1)
    sbs = int(U.mean(0).argmax())
    regret, regret_sbs = (best - U[i, pick]).mean(), (best - U[:, sbs]).mean()
    A = U.shape[1]
    rho = 1 - 6 * ((avg_rank(U) - avg_rank(Up)) ** 2).sum(1) / (A * (A * A - 1))
    return dict(regret=regret, sbs=sbs, regret_sbs=regret_sbs,
                gap=1 - regret / regret_sbs if regret_sbs > 0 else np.nan,
                top1=(U[i, pick] >= best - tol).mean(), rho=rho.mean(), pick=pick)


def probe_first(uv, V, Up, P, cfg, u_fail, S=None, second=None):
    """What the probe rule gets on each mesh: run the probe; keep it if it was valid and ranks first
    in Up; else run `second` (by default the best-ranked other method) and keep the better of the
    two valid runs, and the fallback's run if neither is valid. uv (mesh, method) is the utility of a
    valid run (NaN where none was), V the fraction of runs that were valid (the two runs taken as
    independent when it is a fraction). Returns the utility, the seconds (None without S), whether
    the probe was kept, and the method run second."""
    n, A = uv.shape
    i, z, fb = np.arange(n), cfg["probe"], cfg["fallback"]
    if second is None:
        keep = (P[:, z] >= 0.5) & (Up.argmax(1) == z)
        second = np.where(np.arange(A) == z, -np.inf, Up).argmax(1)
    else:
        keep, second = np.zeros(n, bool), np.full(n, second)
    u = np.where(V > 0, uv, 0.0)
    vz, vx, a, b = V[:, z], V[i, second], u[:, z], u[i, second]
    # Nothing valid kept: the fallback runs, unless it is the second method, which just failed.
    # (A kept probe that is valid on only some copies of the mesh falls back on the others.)
    f_keep = np.full(n, u_fail) if fb < 0 else V[:, fb] * u[:, fb] + (1 - V[:, fb]) * u_fail
    f = np.where(second == fb, u_fail, f_keep)
    ran = vz * vx * np.maximum(a, b) + vz * (1 - vx) * a + (1 - vz) * vx * b + (1 - vz) * (1 - vx) * f
    util = np.where(keep, vz * a + (1 - vz) * f_keep, ran)
    secs = None
    if S is not None:
        s_fb = S[:, fb] if fb >= 0 else np.zeros(n)
        rerun = (1 - vz) * (1 - vx) * np.where(second != fb, s_fb, 0.0)
        secs = S[:, z] + np.where(keep, (1 - vz) * s_fb, S[i, second] + rerun)
    return util, secs, keep, second


def scores(T, V, P, Tp, cfg, S=None):
    """Selection quality of out-of-fold predictions (P, Tp) against the truth (T, V), on meshes
    whose every method has a known validity; S, seconds per run, adds the probe rule's cost."""
    ok = ~np.isnan(V).any(1)
    T, V, P, Tp = T[ok], V[ok], P[ok], Tp[ok]
    S = S[ok] if S is not None else None
    mu, sd = np.nanmean(T, axis=(0, 1)), np.nanstd(T, axis=(0, 1)) + 1e-12
    fb = cfg["fallback"] if cfg["fallback"] >= 0 else None
    z = cfg.get("probe", -1)
    U, Up = utilities(V, T, cfg, mu, sd, None), utilities(P, Tp, cfg, mu, sd, None)
    if z >= 0:      # as Selector.forward: a probe that was not valid is no longer an option
        Up[:, z] = np.where(P[:, z] < 0.5, -np.inf, Up[:, z])
    plain = selection(U, Up)
    s = {"regret": plain["regret"], "gap closed": plain["gap"], "top-1 (ties count)": plain["top1"],
         "rank correlation": plain["rho"]}
    info = dict(n=int(ok.sum()), plain=plain)
    if fb is not None:
        Uf, Upf = utilities(V, T, cfg, mu, sd, fb), utilities(P, Tp, cfg, mu, sd, fb)
        if z >= 0:
            Upf[:, z] = np.where(P[:, z] < 0.5, -np.inf, Upf[:, z])
        with_fb = selection(Uf, Upf)
        s |= {"regret, fallback": with_fb["regret"], "gap closed, fallback": with_fb["gap"]}
        info["fallback"] = with_fb
    if z >= 0:
        k = np.array(cfg["sign"]) * np.array(cfg["weights"]) / sd
        uv = ((T - mu) * k).sum(-1)
        u_fail = float(((np.array(cfg["t_fail"]) - mu) * k).sum())
        # The oracle is the best single method, as for the rows above; the fallback's semantics.
        best = utilities(V, T, cfg, mu, sd, fb).max(1)
        util, secs, keep, second = probe_first(uv, V, Upf if fb is not None else Up, P, cfg, u_fail, S)
        # The rule with no model: the probe and one fixed method, the better kept.
        fixed = [(np.mean(best - probe_first(uv, V, None, P, cfg, u_fail, S, x)[0]), x)
                 for x in range(V.shape[1]) if x != z]
        r_fixed, x_fixed = min(fixed)
        regret = np.mean(best - util)
        name = cfg["algorithms"][z]
        s |= {f"regret, {name} first": regret, f"gap closed, {name} first": 1 - regret / r_fixed if r_fixed > 0 else np.nan,
              f"keeps {name}": keep.mean()}
        if secs is not None:
            s[f"seconds per mesh, {name} first"] = secs.mean()
        pair_secs = probe_first(uv, V, None, P, cfg, u_fail, S, x_fixed)[1]
        info["probe"] = dict(keep=keep, second=second, fixed=x_fixed, regret_fixed=r_fixed,
                             secs_fixed=None if pair_secs is None else pair_secs.mean(),
                             secs_all=None if S is None else S.sum(1).mean())
    # Validity and quality as predictions: the probe's are its result, not predicted, so left out.
    pred = [a for a in range(V.shape[1]) if a != z]
    fails = V < 0.5
    s["validity AUC"] = np.nanmean([auc(1 - P[:, a], fails[:, a]) for a in pred])
    T, Tp = T[:, pred], Tp[:, pred]
    for j, name in enumerate(cfg["metrics"]):
        t, p = T[..., j], Tp[..., j]
        m = ~np.isnan(t)
        with warnings.catch_warnings():   # a mesh no method was valid on has no mean, and no weight
            warnings.simplefilter("ignore", RuntimeWarning)
            tc = np.where(m, t - np.nanmean(t, 1, keepdims=True), 0.0)
        pc = np.where(m, p - np.nansum(np.where(m, p, 0), 1, keepdims=True) / np.maximum(m.sum(1, keepdims=True), 1), 0.0)
        s[f"R2 {name}, between methods"] = 1 - ((tc - pc) ** 2).sum() / (tc ** 2).sum()
    return s, info


# ----------------------------------------------------------------------------------------- main
def make_config(algs, metrics, input_cols, kinds, hidden):
    # A probe column is log1p whatever its training range: a failed probe fills in a 0, which a log
    # would put at -69 if every probe in the training data happened to be valid.
    kinds = [LOG1P if str(c).startswith("probe_") else k for c, k in zip(input_cols, kinds)]
    cfg = dict(algorithms=algs, metrics=metrics, input_cols=input_cols, input_kinds=kinds, hidden=hidden,
               dropout=DROPOUT, n_members=N_SEEDS,
               sign=[1.0 if METRICS[m][0] == "higher" else -1.0 for m in metrics],
               is_log=[bool(METRICS[m][1]) for m in metrics],
               weights=[float(METRICS[m][2]) for m in metrics],
               fallback=algs.index(FALLBACK) if FALLBACK in algs else -1,
               fill=fillable(input_cols), probe=-1, probe_index=[0] * (1 + len(metrics)))
    cfg["t_fail"] = failed_values(cfg)
    if FALLBACK is not None and FALLBACK not in algs:
        print(f"note: FALLBACK '{FALLBACK}' is not one of the methods; a failure is final")
    need = [f"probe_{PROBE}_valid"] + [f"probe_{PROBE}_{m}" for m in metrics]
    if PROBE in algs and all(c in input_cols for c in need):
        cfg.update(probe=algs.index(PROBE), probe_index=[input_cols.index(c) for c in need])
    return cfg


def train(args):
    D, input_cols = load_data(args.data)
    metrics = list(METRICS)
    algs, Y, V = load_targets(D, metrics)
    keep = ~np.isnan(V).all(1)
    if not keep.all():
        print(f"note: {(~keep).sum()} meshes have no results and are skipped")
    D, Y, V = D[keep].reset_index(drop=True), Y[keep], V[keep]
    X = D[input_cols].to_numpy(np.float32)
    X = np.where(np.array(fillable(input_cols)) & np.isnan(X), np.float32(0), X)   # as Selector does
    cfg = make_config(algs, metrics, input_cols, input_kinds(X), HIDDEN)
    T = transform(Y, cfg)
    S = load_seconds(D, algs)
    part_groups = (D["group"] if "group" in D.columns else D["mesh_id"]).astype(str).to_numpy()
    shapes = shape_groups(D)
    n_parts, n_shapes = len(np.unique(part_groups)), len(np.unique(shapes))
    print(f"{len(D)} meshes, {n_shapes} distinct shapes, in {n_parts} groups | {len(input_cols)} inputs "
          f"({', '.join(g for g in INPUTS if input_columns(D, (g,)))}) -> {len(algs)} methods x (valid + "
          f"{len(metrics)} metrics) | methods: {', '.join(algs)} | valid runs: "
          + ", ".join(f"{a} {np.nanmean(V[:, i]):.0%}" for i, a in enumerate(algs))
          + (f" | probe: {PROBE}" if cfg["probe"] >= 0 else ""))
    plot = None
    if args.plot != "off":
        try:
            plot = LossPlot(args.plot_dir, live=args.plot == "live")
        except ImportError:
            print("note: matplotlib is not installed; training without loss plots")

    protocols = {"part": part_groups, "shape": shapes}
    if args.cv != "both":
        protocols = {args.cv: protocols[args.cv]}
    if args.no_cv:
        protocols = {}
    for pname, groups in protocols.items():
        if len(np.unique(groups)) < 2 * N_FOLDS:
            print(f"too few {pname} groups to cross-validate")
            continue
        models = {"linear": []} | ({f"mlp {HIDDEN}": HIDDEN} if HIDDEN else {})
        table = {}
        for name, hidden in models.items():
            t0 = time.time()
            if plot:
                plot.begin_stage(f"cv {pname} {name}")
            P, Tp = cv_predictions({**cfg, "hidden": hidden}, X, T, V, groups, plot)
            table[name], info = scores(T, V, P, Tp, cfg, S)
            if plot:
                plot.end_stage()
            print(f"  cross-validated {name}, {pname} folds, in {time.time() - t0:.0f}s")
        what = ("groups (parts)" if pname == "part" else
                "shapes: meshes with identical boundary features share a fold, so each is tested on shapes it "
                "did not train on")
        print(f"\nOut-of-fold scores, {N_FOLDS}-fold CV by {what}; {N_SEEDS}-seed ensembles, {info['n']} meshes:")
        print(f"  {'':34s}" + "".join(f"{n:>14s}" for n in table))
        for key in table[name]:
            pct = key.startswith(("gap closed", "top-1", "keeps"))
            fmt = "{:>14.0%}" if pct else "{:>14.1f}" if key.startswith("seconds") else "{:>14.3f}"
            print(f"  {key:34s}" + "".join(fmt.format(t[key]) for t in table.values()))
        p = info["plain"]
        print(f"  Regret = utility lost vs. the best method per mesh (0 = oracle), with a failure counting as "
              f"METRICS' failed values. Always picking the single best method ({algs[p['sbs']]}) gives "
              f"{p['regret_sbs']:.3f}; 'gap closed' is the share of that regret removed.")
        if "fallback" in info:
            fb = info["fallback"]
            print(f"  With a failed run replaced by {FALLBACK}'s: always running {algs[fb['sbs']]} (then "
                  f"{FALLBACK} if it fails) gives {fb['regret_sbs']:.3f}, which is the baseline a selector has "
                  f"to beat when failures are detected.")
        if "probe" in info:
            pr = info["probe"]
            print(f"  {PROBE.capitalize()} first = the way to use this model: run {PROBE}; keep it if it was valid and "
                  f"ranks first; else run the best-ranked other method and keep the better of the two"
                  + (f" ({FALLBACK} if neither is valid)" if "fallback" in info else "")
                  + f". With no model, always running {PROBE} + {algs[pr['fixed']]} and keeping the better gives "
                  f"{pr['regret_fixed']:.3f}"
                  + (f" at {pr['secs_fixed']:.1f} s/mesh" if pr["secs_fixed"] is not None else "")
                  + ", the baseline for 'gap closed'"
                  + (f"; all {len(algs)} methods take {pr['secs_all']:.1f} s/mesh." if pr["secs_all"] is not None
                     else "."))
        print("  R2 between methods = variance explained of each metric's differences between the valid "
              "methods on a mesh, which is what the ranking uses"
              + (f"; it and the AUC leave out {PROBE}, whose result is an input." if "probe" in info else "."))
        print(f"  Picked by {name}: " + ", ".join(f"{a} {c:.0%}" for a, c in
                                                zip(algs, np.bincount(p["pick"], minlength=len(algs)) / info["n"])))
        if "probe" in info:
            pr = info["probe"]
            runs = np.bincount(pr["second"][~pr["keep"]], minlength=len(algs)) / info["n"]
            print(f"  {PROBE.capitalize()} first, {name}: keeps {PROBE} {pr['keep'].mean():.0%}, else runs "
                  + ", ".join(f"{a} {c:.0%}" for a, c in zip(algs, runs) if a != PROBE))
        # On every regret row: one alone moves by ~0.01 with the fold assignment, which decides nothing.
        if len(table) == 2 and all(table["linear"][k] <= table[name][k] for k in table[name] if k.startswith("regret")):
            print("  The linear model is as good as the MLP: prefer it (HIDDEN = []) or try LightGBM.")
        print()

    t0 = time.time()
    if plot:
        plot.begin_stage("final")
    sel = fit_selector(cfg, X, T, V, part_groups, seed=999, plot=plot)
    if plot:
        plot.end_stage()
    torch.save({"config": cfg, "state": sel.state_dict()}, args.model)
    print(f"Trained on all {len(D)} meshes in {time.time() - t0:.0f}s -> {args.model}")
    if args.onnx:
        export_onnx(sel, cfg, args.onnx)


def export_onnx(sel, cfg, path):
    x = torch.ones(2, len(cfg["input_cols"]))
    names = dict(input_names=["x"], output_names=["metrics", "valid", "utility"])
    logging.getLogger("torch.onnx").setLevel(logging.ERROR)
    try:      # torch.export-based exporter (PyTorch >= 2.5; needs the onnxscript package)
        torch.onnx.export(sel, (x,), path, dynamo=True, verbose=False,
                          dynamic_shapes=({0: torch.export.Dim("n")},), **names)
    except Exception:
        legacy = dict(dynamic_axes={"x": {0: "n"}, "metrics": {0: "n"}, "valid": {0: "n"}, "utility": {0: "n"}},
                      **names)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                torch.onnx.export(sel, (x,), path, dynamo=False, **legacy)
            except TypeError:        # PyTorch < 2.5 has only the legacy exporter
                torch.onnx.export(sel, (x,), path, **legacy)
    A, M = len(cfg["algorithms"]), len(cfg["metrics"])
    z = cfg.get("probe", -1)
    print(f"ONNX -> {path}\n  input  x: float32 [n, {x.shape[1]}] = {cfg['input_cols'][0]}..{cfg['input_cols'][-1]}"
          + (f" ({cfg['algorithms'][z]}'s metrics and block count may be NaN)" if z >= 0 else "") + "\n"
          f"  output metrics: [n, {A}, {M}] (methods {cfg['algorithms']}, metrics {cfg['metrics']}), of a valid run\n"
          f"  output valid:   [n, {A}], the probability of a valid decomposition\n"
          f"  output utility: [n, {A}], higher = better"
          + (f"; keep {cfg['algorithms'][z]} if it is first, else run the first and keep the better" if z >= 0 else ""))


def predict(args):
    ck = torch.load(args.model, map_location="cpu", weights_only=True)
    cfg = ck["config"]
    if "input_cols" not in cfg:
        sys.exit(f"{args.model} is from before the validity head and the feature inputs; retrain it")
    sel = Selector(cfg)
    sel.load_state_dict(ck["state"])
    sel.eval()
    if set(METRICS) == set(cfg["metrics"]):   # use the current utility weights and fallback
        sel.w.copy_(torch.tensor([float(METRICS[m][2]) for m in cfg["metrics"]]))
        sel.fallback.fill_(cfg["algorithms"].index(FALLBACK) if FALLBACK in cfg["algorithms"] else -1)
    D, input_cols = load_data(args.data, cfg["input_cols"])
    X = torch.tensor(D[input_cols].to_numpy(np.float32))    # a missing probe metric: the model fills it
    with torch.no_grad():
        Y, P, U = sel(X)
    order = U.argsort(1, descending=True).numpy()
    algs, z = cfg["algorithms"], cfg.get("probe", -1)
    # The probe rule: keep the probe when it ranks first (only a valid one can), else run the
    # best-ranked other method. None = nothing more to run.
    run = [None if z < 0 or o[0] == z else int(o[0]) for o in order]
    rows = [{"mesh_id": mid, "rank": r + 1, "method": algs[a], "utility": float(U[i, a]), "p_valid": float(P[i, a]),
             **({"action": "keep" if a == z and run[i] is None else "run" if a == run[i] else ""} if z >= 0 else {}),
             **{m: float(Y[i, a, j]) for j, m in enumerate(cfg["metrics"])}}
            for i, mid in enumerate(D["mesh_id"]) for r, a in enumerate(order[i])]
    pd.DataFrame(rows).to_csv(args.out, index=False)
    for i, mid in enumerate(D["mesh_id"][:10]):
        ranking = " > ".join(f"{algs[a]} ({P[i, a]:.0%} valid)" for a in order[i])
        if z < 0:
            print(f"{mid}: {ranking}")
            continue
        ok = P[i, z] >= 0.5
        then = "keep it" if run[i] is None else f"run {algs[run[i]]}" + (", keep the better" if ok else "")
        print(f"{mid}: {algs[z]} {'valid' if ok else 'not valid'} -> {then}   [{ranking}]")
    weights = ", ".join(f"{m} {w:g}" for m, w in zip(cfg["metrics"], sel.w.tolist()))
    fb = int(sel.fallback)
    print(("...\n" if len(D) > 10 else "") + f"Utility weights: {weights}; "
          + (f"a failed run replaced by {algs[fb]}'s" if fb >= 0 else "a failure is final")
          + (f"; {algs[z]} has run, and 'action' says whether to keep it or which method to run next "
             f"(then keep the better of the two)" if z >= 0 else "")
          + f". Rankings, validity and predicted metrics for {len(D)} meshes -> {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="dataset.csv from build_dataset.py (CSV or Parquet)")
    ap.add_argument("--model", default="selector.pt", help="model file to save, or to load with --predict")
    ap.add_argument("--predict", action="store_true", help="rank the methods for the meshes in --data")
    ap.add_argument("--out", default="predictions.csv", help="output file for --predict")
    ap.add_argument("--no-cv", action="store_true", help="skip cross-validation; just train and save")
    ap.add_argument("--cv", choices=["both", "part", "shape"], default="both",
                    help="cross-validate by group (part), by identical boundary features (shape), or both")
    ap.add_argument("--plot", choices=["live", "save", "off"], default="live",
                    help="loss curves while training: a live window + PNGs, PNGs only, or none")
    ap.add_argument("--plot-dir", default="loss_plots", help="where the loss_<stage>.png files go")
    ap.add_argument("--onnx", help="also export the trained model to this ONNX file")
    args = ap.parse_args()
    predict(args) if args.predict else train(args)


if __name__ == "__main__":
    main()
