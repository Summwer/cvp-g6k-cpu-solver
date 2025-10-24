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

vector<vector<pii>> make_graph_cvp(double a, int grid, double c=1);
double dijkstra( vector<vector<pii>> G, int s, int t );
int predict_slicer_max_sample_amount(double a, double len_bound, int n);