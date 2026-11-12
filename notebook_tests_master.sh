#!/bin/bash

# Run by MuPT's CI (it clones this repo and sources this script), so this list
# keeps testing the same notebooks as before the reorganization.
# workflows/ is not listed: it needs flowerMD + HOOMD, which MuPT's CI env lacks.
example_dirs=(
    'representation'
    'interfaces'
    'legacy'
)
batch_status="passing"

for dirname in "${example_dirs[@]}"; do
    pushd "./$dirname"
    bash "./notebook_tests.sh" || batch_status="failing"
    popd
done

if [ $batch_status = "failing" ]; then
    exit 1 # bubble up failure to CI if any of batches fail
fi