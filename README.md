# cvp-g6k-cpu-solver


I implement the slicer algorithm in g6k kernel in cpp, and it allows to call it in python. I also implemented the progressive slicer algorithm in python file named "pro_randslicer.py". 



first compile command:
```
PYTHON=python2 ./bootstrap.sh
#or
PYTHON=python3 ./bootstrap.sh
```


recompile command:
```
python setup.py clean
python setup.py build_ext --inplace 
```


### Compute Expected Sample times for Randomized Slicer with d4f

```
cd kernel
g++ -o cm compute_acvp_sample_estimation.cpp  slicer_sample_estimation.cpp
./cm
```



### Test for Slicer

In `slicer_test.py`, update the parameter setting of `max_sample_times` and `len_bound` in this file to help other developers adjust them. 

Besides, we update the test result of slicer while setting `max_sample_times = 140`.


Test time cost of one iterative slicer tour and whether it can find a close vector whose distance is sqrt(4/3)gh within target vector t.
```
python -u slicer_test.py --threads 1 --max-sample-times 1 --approx-factor 1.1547 --max-dim 95 |& tee test_results/slicer_tests/iterative-slicer-prebkz-bdgl2-threads=1.log  2>&1
```

Test time cost of one randomized iterative slicer tour and whether it can find a close vector whose distance is gh within target vector t.
```
python -u slicer_test.py --approx-factor 1.1  --max-dim 90 |& tee test_results/slicer_tests/random-slicer-prebkz-bdgl2-threads=1.log  2>&1
```


Test time cost of one randomized iterative slicer tour and whether it can find a close vector whose distance is gh within target vector t using 20 threads.
```
python -u slicer_test.py --threads 20 --approx-factor 1.1 --max-dim 90 |& tee test_results/slicer_tests/random-slicer-prebkz-bdgl2-threads=20.log  2>&1
```

#### Test slicer for batch-CVP with small approximate factor

Test time cost of one randomized iterative slicer tour to solve a batch-CVP problem with batch size 100, whose distance is gh within target vector t using 20 threads.
```
python -u slicer_test.py --threads 20 --approx-factor 1.1 --batch-size 100 --max-dim 90  |& tee test_results/slicer_tests/random-slicer-prebkz-bdgl2-batch-cvpthreads=20.log  2>&1
```



#### Test the d4f techinque for randomized slicer


Test the correctness of the d4f funtion in [Leo18] in slicer. 
```
python -u slicer_test.py  --approx-factor 1.1547 --batch-size 1 --max-dim 90 --consider-d4f True --max-sample-times 1000 | tee test_results/slicer_tests/d4f-test/test-cvp-d4f-function.log 2>&1
```


Test the correctness of different d4f values in slicer. 
./test-d4f.sh
```







### Test for Colattice

#### Test Colattice with Trivial Strategy

Test time cost of one colattice tour with original colattice strategy. Find a close vector whose distance is 2.2*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.2 --min-dim 150 --max-dim 155 |& tee test_results/colattice_tests/trival-colattice-threads=1-22.log  2>&1
```


Test time cost of one colattice tour with original colattice strategy. Find a close vector whose distance is 2.5*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.5 --min-dim 160 --max-dim 175 |& tee test_results/colattice_tests/trival-colattice-threads=1-25.log  2>&1
```

Test time cost of one approximate colattice tour with original colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.2*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.2 --min-dim 150 --max-dim 155 --len-bound 1.333 |& tee test_results/colattice_tests/trival-approx-colattice-threads=1-22-1333.log  2>&1
```


Test time cost of one approximate colattice tour with original colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.5*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.5 --min-dim 160 --max-dim 175 --len-bound 1.333 |& tee test_results/colattice_tests/trival-approx-colattice-threads=1-25-1333.log  2>&1
```




#### Test Colattice with Optimized Strategy


Test time cost of one colattice tour with optimized colattice strategy and whether it can find a close vector. Find a close vector whose distance is 2.2*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.2 --min-dim 150 --max-dim 155 --optimize-strategy True |& tee test_results/colattice_tests/op-colattice-threads=1-22.log  2>&1
```

Test time cost of one colattice tour with optimized colattice strategy and whether it can find a close vector. Find a close vector whose distance is 2.5*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.5 --min-dim 160 --max-dim 175 --optimize-strategy True |& tee test_results/colattice_tests/op-colattice-threads=1-25.log  2>&1
```


Test time cost of one approximate colattice tour with optimized colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.2*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.2 --min-dim 150 --max-dim 155 --len-bound 1.333 --optimize-strategy True |& tee test_results/colattice_tests/op-approx-colattice-threads=1-22-1333.log  2>&1
```

Test time cost of one approximate colattice tour with optimized colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.5*gh within target vector t. 

```
python -u colattice_test.py --approx-factor 2.5 --min-dim 160 --max-dim 175 --len-bound 1.333 --optimize-strategy True |& tee test_results/colattice_tests/op-approx-colattice-threads=1-25-1333.log  2>&1
```


#### Test for batch-CVP



Test time cost of one approximate colattice tour with optimized colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.2*gh within target vectors with batch size 100 using 20 threads.

```
python -u colattice_test.py --approx-factor 2.2 --min-dim 150 --max-dim 155 --len-bound 1.333 --batch-size 100 --threads 20 --optimize-strategy True |& tee test_results/colattice_tests/op-approx-colattice-threads=20-22-1333-batch-cvp.log  2>&1
```

Test time cost of one approximate colattice tour with optimized colattice strategy, gamma1 = sqrt(4/3). Find a close vector whose distance is 2.5*gh within target vectors with batch size 100 using 20 threads.

```
python -u colattice_test.py --approx-factor 2.5 --min-dim 160 --max-dim 175 --len-bound 1.333 --batch-size 100 --threads 20 --optimize-strategy True |& tee test_results/colattice_tests/op-approx-colattice-threads=20-25-1333-batch-cvp.log  2>&1
```