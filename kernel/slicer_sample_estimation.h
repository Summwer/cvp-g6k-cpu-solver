#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <math.h>
#include <fstream>
#include <omp.h>

using namespace std;

typedef pair<double, int> pii;
typedef pair<double, pair<int,int>> qpii;


int predict_slicer_max_sample_amount(double a, double len_bound, int n);