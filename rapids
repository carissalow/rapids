#!/usr/bin/env python
from sys import argv
import subprocess
subprocess.run(" ".join(["snakemake", "-R", "`snakemake --list-params-changes`"] + argv[1:]), shell=True)