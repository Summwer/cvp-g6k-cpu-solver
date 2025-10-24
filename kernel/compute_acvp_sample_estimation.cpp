

#include "slicer_sample_estimation.h"
#include <iostream>

int main(){
    double len_bound = 1.;
    double saturation_ratio = 0.5;
    double db_size_base = pow((4./3.), 0.5);
    int actual_n;

    for(int len_bound_step = 0; len_bound_step < 5; len_bound_step++){
        len_bound = 1.05 + len_bound_step * 0.01;
        std::cout << "len_bound = " << len_bound << std::endl;
        int grid = 5000;
        auto G = make_graph_cvp( sqrt(4/3.), grid, len_bound );
        double plog = dijkstra(G, grid, 0);
        std::cout << "plog = " << plog << std::endl;
        std::cout<<"{";
        for(int n = 60; n <= 85; n +=5 ){
            if(n==60 or n == 65)
                actual_n = n - 13;
            else if(n==70 or n==75 or n== 80)
                actual_n = n - 14;
            else if(n==85)
                actual_n = n - 15;
            double a = pow(.5 * saturation_ratio, 1./actual_n ) + db_size_base;
            double pn = predict_slicer_max_sample_amount(a, pow(len_bound,2), actual_n );
            std::cout << n << ":" << pn << ", ";
        }
        std::cout<<"}"<<std::endl;
    }
    
}
