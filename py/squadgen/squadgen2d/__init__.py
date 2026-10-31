"""
squadgen2d -- SQuadGen (docs/squadgen_2d.md) for planar CrossGen domains.

The pipeline, one command per stage (python -m squadgen2d <command>, or `squadgen2d <command>` once installed with
`pip install -e py/squadgen`):

    cdf-data     CrossGen layouts -> CDF NPZs (cdf_data.py; the only stage that needs the C++ crossgen module)
    cdf-pack     CDF NPZs -> presaved fields at every training point (cdf_pack.py)
    geom-ae      OBJ domains -> Geometry Autoencoder (geom_ae.py)
    sq-vae       field pack + frozen Geom-AE -> SQ-VAE (sq_vae.py)
    sq-diffuse   field pack + frozen SQ-VAE -> SQ-Diffuse (sq_diffuse.py)
    roundtrip    SQ-VAE -> QuadExtraction on true vs decoded fields (roundtrip.py)

The library modules under them:

    paths        the repository, py/squadgen/{data,runs}, and repository-relative checkpoint paths
    training     device, mixed precision, optimiser step, DataLoader, schedule and log, shared by the trainers
    geometry     numpy geometry: OBJ, sampling, FPS, truncated SDF, frames, Deformation (augmentation)
    domains      the Geom-AE's domains, families and split, and the FPS token-set cache
    layout       a CDF NPZ's split mesh and its fields (Appendix B)
    layers       the network modules the Geom-AE and SQ-VAE share

Importing the package imports none of them, so `import squadgen2d.geometry` does not pull in torch. Models load with
geom_ae.load_geom_ae(), sq_vae.load_sq_vae() and sq_diffuse.load_sq_diffuse().
"""

__version__ = "0.2.0"
