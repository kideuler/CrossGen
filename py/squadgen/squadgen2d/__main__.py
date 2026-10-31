"""python -m squadgen2d <command> [options]  (or `squadgen2d <command>` once installed): see COMMANDS."""

import importlib
import sys

# command -> (module, what it does). Each module's main(argv) takes the command's options; `<command> -h` lists them.
COMMANDS = {
    "cdf-data": ("cdf_data", "make CDF NPZs of CrossGen layouts (needs the crossgen module)"),
    "cdf-pack": ("cdf_pack", "presave a cdf-data directory's fields for training"),
    "geom-ae": ("geom_ae", "train / evaluate the Geometry Autoencoder"),
    "sq-vae": ("sq_vae", "train / evaluate the SQ-VAE on a field pack"),
    "sq-diffuse": ("sq_diffuse", "train / evaluate SQ-Diffuse"),
    "roundtrip": ("roundtrip", "SQ-VAE round trip through QuadExtraction"),
}


def usage():
    width = max(map(len, COMMANDS))
    return ("usage: python -m squadgen2d <command> [options]\n\ncommands, in pipeline order:\n"
            + "\n".join(f"  {c:{width}s}  {what}" for c, (_, what) in COMMANDS.items())
            + "\n\n`python -m squadgen2d <command> -h` lists a command's options.")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(usage())
        return 0
    if argv[0] not in COMMANDS:
        print(f"unknown command {argv[0]!r}\n\n{usage()}", file=sys.stderr)
        return 2
    module = importlib.import_module(f"{__package__ or 'squadgen2d'}.{COMMANDS[argv[0]][0]}")
    return module.main(argv[1:]) or 0


if __name__ == "__main__":
    raise SystemExit(main())
