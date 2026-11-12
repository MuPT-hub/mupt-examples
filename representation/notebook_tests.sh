#!/bin/bash

py.test --nbval-lax --nbval-current-env -vv \
    01_hierarchy_basics.ipynb \
    02_hierarchy_advanced.ipynb \
    03_shapes.ipynb \
    04_connector_alignment.ipynb
