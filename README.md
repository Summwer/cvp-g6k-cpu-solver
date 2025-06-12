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

In `slicer_test.py`, update the parameter setting of `max_sample_times` and `len_bound` in this file to help other developers adjust them. 

Besides, we update the test result of slicer while setting `max_sample_times = 140` and use the cvp instance generation method from https://github.com/RandomizedSlicer/RandomizedSlicerG6K. Under the same benchmark, our running time is shorter than it.  