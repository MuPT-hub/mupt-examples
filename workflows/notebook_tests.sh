#!/bin/bash
# Workflow notebooks need flowerMD + HOOMD (conda-envs/workflows.yml), which
# MuPT's CI environment does not have yet, so notebook_tests_master.sh does not
# run this batch. Run it by hand from this folder in the workflows environment.

py.test --nbval-lax --nbval-current-env -vv \
    01_pe_all_atom_dpd.ipynb
