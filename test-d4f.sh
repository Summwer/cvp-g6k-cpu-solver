#!/bin/bash

# Loop through d4f values from 1 to 10
for d4f in {1..20}
do
    python -u slicer_test.py --approx-factor 1.1547 --batch-size 1 --max-dim 87 --consider-d4f True --d4f $d4f --max-sample-times 1000 | tee test_results/slicer_tests/d4f-test/test-cvp-d4f-$d4f.log 2>&1 &
done

# Wait for all background jobs to finish
wait