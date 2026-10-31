"""
What the three trainers share: the device, mixed precision, the optimiser step, the DataLoader, the learning-rate
schedule and the run log.

Devices. pick_device() takes --device or the first of CUDA, Apple MPS and the CPU. On an NVIDIA GPU the forward pass
runs under bf16 autocast (--amp auto; fp16 with loss scaling for a GPU without bf16), as the release trains, and
float32 matmuls may use TF32; on MPS autocast was slower, so it stays fp32 there. The weights are fp32 either way and
checkpoints hold them on the CPU, so a model trained on a cloud GPU loads and runs on a Mac as it is.

The release trainers need CUDA and torchrun (docs/squadgen_2d.md Sec. 3.4, "Hardware"), so each trainer here is the
plain loop around its loss that the guide suggests.
"""

import contextlib
import math
import os
import time

import torch


def pick_device(pref=None):
    """The device to run on: `pref` when given ('cuda', 'cuda:1', 'mps', 'cpu'), else the first of CUDA, Apple MPS
    and the CPU that this machine has. On CUDA it also lets float32 matmuls use TF32 (Ampere and later), several
    times faster at a precision these models do not notice; MPS and the CPU have no such mode."""
    if pref and pref != "auto":
        dev = pref
        if dev.startswith("cuda") and not torch.cuda.is_available():
            raise SystemExit(f"--device {dev}: this machine has no CUDA device")
        if dev == "mps" and not torch.backends.mps.is_available():
            raise SystemExit("--device mps: this machine has no MPS device")
    else:
        dev = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
    if dev.startswith("cuda"):
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
    return dev


def device_type(dev):
    return dev.split(":")[0]


def amp_dtype(dev, mode):
    """The autocast dtype of the training forward pass, or None for fp32. 'auto' is bf16 on a CUDA device that has
    it (Ampere and later: A100, H100, RTX 30xx/40xx/50xx), as the release trains (--mix_precision_dtype=bf16), and
    fp32 elsewhere: on MPS bf16 autocast measured slower, and a T4 or V100 has no fast bf16, for which 'fp16' (with
    loss scaling) is the choice. Only the arithmetic changes: the weights stay fp32, so whatever trained under
    autocast runs in fp32 on any device. Evaluation always runs in fp32."""
    if mode == "none":
        return None
    if mode == "auto":
        return torch.bfloat16 if dev.startswith("cuda") and torch.cuda.is_bf16_supported() else None
    return {"bf16": torch.bfloat16, "fp16": torch.float16}[mode]


def autocast(dev, dtype):
    return contextlib.nullcontext() if dtype is None else torch.autocast(device_type(dev), dtype=dtype)


def grad_scaler(dev, dtype):
    """fp16 needs its loss scaled so that small gradients do not underflow; bf16 and fp32 do not, and a disabled
    scaler passes everything through."""
    return torch.amp.GradScaler(device_type(dev), enabled=dtype == torch.float16)


def train_step(model, opt, scaler, loss, clip):
    """Backward, gradient clipping at `clip`, optimiser step; returns the gradient norm (inf on an fp16 step the
    scaler skipped)."""
    opt.zero_grad(set_to_none=True)
    scaler.scale(loss).backward()
    scaler.unscale_(opt)
    gnorm = torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
    scaler.step(opt)
    scaler.update()
    return gnorm


def lr_at(step, cfg):
    """Linear warm-up, then half-cycle cosine to min_lr (the release's lr_sched, per step)."""
    warm = max(1, int(cfg["warmup_frac"] * cfg["steps"]))
    if step < warm:
        return cfg["lr"] * (step + 1) / warm
    t = (step - warm) / max(1, cfg["steps"] - warm)
    return cfg["min_lr"] + (cfg["lr"] - cfg["min_lr"]) * 0.5 * (1 + math.cos(math.pi * t))


def accumulate(run, terms):
    """Sum the step's logged terms into `run` (key -> (sum, steps)) until the next log line. A gradient norm that is
    not finite is counted instead of summed: fp16 steps the scaler skips have one, and averaging it in would hide
    every other step's. Any other term is summed as it is, so a loss that has gone NaN shows."""
    for k, v in terms.items():
        x = float(v.detach()) if torch.is_tensor(v) else float(v)
        if k == "gnorm" and not math.isfinite(x):
            run["nonfinite"] = run.get("nonfinite", 0) + 1
            continue
        s, n = run.get(k, (0.0, 0))
        run[k] = (s + x, n + 1)


def averages(run):
    out = "  ".join(f"{k} {s / n:.4f}" for k, v in run.items() if k != "nonfinite" for s, n in [v])
    return out + (f"  [{run['nonfinite']} step(s) with a non-finite gradient norm]" if run.get("nonfinite") else "")


def default_workers(dev):
    """DataLoader workers when --workers is not given. Each builds samples in numpy and holds its own torch: two keep
    an M2 fed at batch 8 within 8 GB, while a GPU takes a step in a fraction of that and needs more."""
    return min(8, max(2, (os.cpu_count() or 2) - 1)) if dev.startswith("cuda") else 2


def make_loader(dataset, cfg):
    w = cfg["workers"]
    return torch.utils.data.DataLoader(dataset, batch_size=cfg["batch"], num_workers=w, persistent_workers=w > 0,
                                       prefetch_factor=4 if w > 0 else None,
                                       pin_memory=cfg["device"].startswith("cuda"))


def to_dev(batch, dev):
    return {k: torch.as_tensor(v).to(dev, non_blocking=True) for k, v in batch.items()}


def memory_note(dev):
    """The device memory in use, for the step log ('' on the CPU)."""
    if dev == "mps":
        return f", mps {torch.mps.driver_allocated_memory() / 2 ** 30:.2f} GiB"
    if dev.startswith("cuda"):
        return f", cuda peak {torch.cuda.max_memory_allocated(dev) / 2 ** 30:.2f} GiB"
    return ""


def empty_cache(dev):
    if dev == "mps":
        torch.mps.empty_cache()
    elif dev.startswith("cuda"):
        torch.cuda.empty_cache()


def add_device_args(ap):
    """The options every trainer takes the same way."""
    ap.add_argument("--workers", type=int, help="DataLoader workers (default: 2 on MPS/CPU, up to 8 on CUDA)")
    ap.add_argument("--device", default="auto", help="cuda, cuda:N, mps or cpu (default: the first there is)")
    ap.add_argument("--amp", choices=("auto", "bf16", "fp16", "none"), default="auto",
                    help="training autocast (auto: bf16 on CUDA that has it, else fp32)")
    ap.add_argument("--log_every", type=int, default=100)
    ap.add_argument("--eval_every", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)


class RunLog:
    """print() that also appends to <out>/log.txt, flushed per line, so that a run's log survives a crash and can be
    followed with tail -f."""

    def __init__(self, out):
        os.makedirs(out, exist_ok=True)
        self.file = open(os.path.join(out, "log.txt"), "a")

    def __call__(self, s):
        print(s, flush=True)
        self.file.write(s + "\n")
        self.file.flush()

    def header(self, cfg, skip=()):
        """The run's first line: the time and every config entry, as each run of a directory appends its own."""
        self(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  "
             + " ".join(f"{k}={v}" for k, v in sorted(cfg.items()) if k not in skip))
