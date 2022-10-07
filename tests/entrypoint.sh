#!/bin/sh -l

Rscript code/run_tests.R
RUNNING_TESTS=Y python3 tests/*.py

