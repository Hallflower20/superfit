#!/usr/bin/env python
"""Run NGSF from a source checkout.

Kept so `python run.py parameters.json` keeps working. The implementation
lives in NGSF/cli.py, which is also installed as the `superfit` command.
"""

import sys

from NGSF.cli import main

if __name__ == "__main__":
    sys.exit(main(prog="run.py"))
