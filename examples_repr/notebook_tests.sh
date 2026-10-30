#!/bin/bash

py.test --nbval-lax --current-env -vv \
    mupt_shapes.ipynb \
    connector_alignment.ipynb \
    mupt_primitives_rdkit.ipynb \
    mol_from_scratch_basic.ipynb \
    mol_from_scratch_advanced.ipynb \
