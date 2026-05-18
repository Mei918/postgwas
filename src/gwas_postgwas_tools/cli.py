from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional

from gwas_postgwas_tools.doctor import run_doctor


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gwas-postgwas",
        description="Run post-GWAS validation and analysis workflows.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate-config", help="Validate a YAML config.")
    validate.add_argument("--config", required=True, help="Path to a YAML config file.")

    run = subparsers.add_parser("run-postgwas", help="Run the post-GWAS workflow.")
    run.add_argument("--config", required=True, help="Path to a YAML config file.")

    subparsers.add_parser("doctor", help="Check server-side CLI and dependency availability.")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)

    if raw_argv and raw_argv[0] == "doctor":
        return run_doctor(raw_argv[1:])

    parser = build_parser()
    args = parser.parse_args(raw_argv)

    try:
        if args.command == "validate-config":
            from gwas_postgwas_tools.config import config_to_dict, load_config

            config = load_config(args.config)
            payload = config_to_dict(config)
            print(f"Config is valid for study: {payload['study']['name']}")
            print(f"eqtl_path\t{config.eqtl.path}")
            print(f"gwas_path\t{config.gwas.path}")
            print(f"batch_mode\t{config.batch_mode}")
            return 0

        if args.command == "run-postgwas":
            from gwas_postgwas_tools.config import load_config
            from gwas_postgwas_tools.workflows.postgwas import run_postgwas_pipeline

            config = load_config(args.config)
            outputs = run_postgwas_pipeline(config)
            for label, path in outputs.items():
                print(f"{label}\t{Path(path).resolve()}")
            return 0

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    parser.print_help(sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
