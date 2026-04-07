#!/usr/bin/env python3
"""
TACO - Telomere-Aware Contig Optimization

Entry point for running TACO as: python3 -m taco [options]
"""
import os
import sys

# Ensure the taco package is importable from the script's parent directory
TACO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if TACO_DIR not in sys.path:
    sys.path.insert(0, TACO_DIR)

from taco.cli import parse_args
from taco.pipeline import PipelineRunner


def main():
    args = parse_args()
    runner = PipelineRunner(args)
    runner.run()


if __name__ == "__main__":
    main()
