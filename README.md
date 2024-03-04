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