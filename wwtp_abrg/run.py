"""Command-line entrypoint."""

from __future__ import annotations

import argparse
import logging

from wwtp_abrg.config import load_config
from wwtp_abrg.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the WWTP ABRG pipeline")
    parser.add_argument("--config", required=True, help="Path to YAML configuration file")
    return parser


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = build_parser()
    args = parser.parse_args()
    config = load_config(args.config)
    run_pipeline(config)


if __name__ == "__main__":
    main()
