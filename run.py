#!/usr/bin/env python
"""Run superfit from a source checkout.

Kept so `python run.py parameters.json` keeps working. The implementation
lives in superfit/cli.py, which is also installed as the `superfit` command.
"""

import sys

from superfit.cli import main

if __name__ == "__main__":
    sys.exit(main(prog="run.py"))
