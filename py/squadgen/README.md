# squadgen2d

SQuadGen (Geom-AE → SQ-VAE → SQ-Diffuse, `docs/squadgen_2d.md`) for planar CrossGen domains, as one Python package.

```
py/squadgen/
  pyproject.toml
  squadgen2d/        the package
  data/              CDF data (cdf-data output, field packs)        -- generated, untracked
  runs/              training runs (checkpoints, logs, plots)       -- generated, untracked
```

## Running

Install it once (editable, so edits take effect without reinstalling), then use the `squadgen2d` command:

```sh
python3 -m pip install -e py/squadgen
squadgen2d                       # lists the commands
squadgen2d sq-vae -h             # a command's options
```

Or, with nothing installed, put the project on the path and use `python -m`:

```sh
PYTHONPATH=py/squadgen python3 -m squadgen2d sq-vae -h
```

`squadgen2d <command>` and `python -m squadgen2d <command>` are the same thing; so is `python -m squadgen2d.sq_vae`.
Defaults are paths in this checkout, so commands can run from any directory.

## The pipeline

| Command | In → out | Needs |
|---|---|---|
| `cdf-data run [DIR...]` | OBJ domains → one CDF NPZ per layout + `manifest.csv` | the C++ `crossgen` module |
| `cdf-pack [DATA]` | a cdf-data directory → its presaved fields, `DATA/pack_N512/` | numpy |
| `geom-ae` | OBJ domains → Geometry Autoencoder | torch |
| `sq-vae` | field pack + frozen Geom-AE → SQ-VAE, `latent_stats.json` | torch |
| `sq-diffuse` | field pack + frozen SQ-VAE → SQ-Diffuse | torch |
| `roundtrip` | SQ-VAE → QuadExtraction on true vs decoded fields | QuadTools |

A full run on the ORACLE layouts:

```sh
squadgen2d cdf-data run --methods oracle --out py/squadgen/data/cdf_oracle     # local: needs crossgen
squadgen2d cdf-pack py/squadgen/data/cdf_oracle                                # 5 s, 1 GB
squadgen2d geom-ae --out py/squadgen/runs/geom_ae_2d
squadgen2d sq-vae --data py/squadgen/data/cdf_oracle --geom_ae py/squadgen/runs/geom_ae_2d/best.pt \
    --out py/squadgen/runs/sq_vae_2d_oracle
squadgen2d sq-diffuse --vae py/squadgen/runs/sq_vae_2d_oracle/best.pt --data py/squadgen/data/cdf_oracle
```

`cdf-data` is the only step that needs CrossGen. Copy its output directory to a GPU box and build the pack there:
packing needs only numpy.

## One layout per domain (`cdf-data run --fallback`)

By default `cdf-data run` keeps every method's layout of a domain. With `--fallback` it keeps one:

```sh
squadgen2d cdf-data run data/objs/2d_data --fallback zipline atlas     # -> py/squadgen/data/cdf_zipline_atlas
```

ZIPLINE runs first. Its layout is kept when the blocks are valid with 100 % coverage and their CDF mesh does not fold.
Otherwise ATLAS runs and is held to the same test. A domain neither passes is skipped: it is printed as it happens,
with each method's reason, and listed in `skipped.txt` of the output directory. `--fallback` with no methods means
`zipline atlas`; any chain of methods works. Each method is its own task, so the run resumes after Ctrl-C like any other,
and a directory holds one kind of run: the command refuses to mix chains, or a chain and a run of every method, in one.

## Presaved fields (`cdf-pack`)

A CDF NPZ holds a layout's split mesh, and the fields have to be evaluated on that mesh. A pack evaluates them once, at
every point training reads: the 5 token sets of N and 4N points and the 20 000 query samples of each layout. It stores
them as float32 arrays (`tok512.npy`, `tok2048.npy`, `query.npy`, `boundary.npy`, indexed by `layouts.csv`).
Training memory-maps them, so a sample is a few slices instead of rebuilding the mesh, which was 11 of a sample's 16 ms.

- Unaugmented samples are bit-identical to the NPZ path, so evaluations are unchanged.
- Augmentation moves the presaved points through the same random deformation as before, with no mesh: positions and
  corners are deformed, gradients turned by the Jacobian, and sdf is measured to the deformed boundary. cdf and dcdf
  do not change under a deformation. The result is within 6e-4 of the mesh version in the radius-1 frame.
- A pack is rebuilt only when the manifest's layouts change (or with `--force`), and the trainers refuse a stale one.
- `--tokens` must be the Geom-AE's N (512 by default).

## Train / validation split

The Geom-AE holds out 10 % of the families (CAD models) and records them in its checkpoint. SQ-VAE validates on
layouts of those families, and records them too; SQ-Diffuse and `roundtrip` use the SQ-VAE's. That way no network in
the chain trained on a validation domain, whatever data the Geom-AE was trained on. Older Geom-AE checkpoints did not
record the families, so they are read from the run's `log.txt`. The manifest's `split` column is only a fallback.

## Importing

```python
from squadgen2d.geom_ae import load_geom_ae, encode_obj
from squadgen2d.sq_vae import load_sq_vae, encode_layout, decode_fields
from squadgen2d.sq_diffuse import load_sq_diffuse, domain_condition, sample_latents
from squadgen2d.cdf_pack import FieldPack, make_sample
from squadgen2d import geometry, layout, paths
```

`import squadgen2d` loads nothing. `geometry`, `domains`, `layout`, `cdf_pack` and `cdf_data` need only numpy.

| Module | What |
|---|---|
| `paths` | repository paths, defaults, repo-relative checkpoint paths |
| `training` | device, mixed precision, optimiser step, DataLoader, lr schedule, run log |
| `geometry` | OBJ, area/FPS sampling, truncated SDF, frames, `Deformation` (augmentation) |
| `domains` | the Geom-AE's domains, families, split, FPS token-set cache |
| `layout` | a CDF NPZ's split mesh and its fields (Appendix B) |
| `layers` | network modules shared by the Geom-AE and SQ-VAE |
| `cdf_data`, `cdf_pack`, `geom_ae`, `sq_vae`, `sq_diffuse`, `roundtrip` | the commands |
