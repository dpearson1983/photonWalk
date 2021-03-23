#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include "photon.hpp"

std::vector<double> stats(std::vector<double> &data) {
    double mean = 0.0, M2 = 0.0;
    int count = 0;
    for (int i = 0; i < data.size(); ++i) {
        count++;
        double delta = data[i] - mean;
        mean += delta/count;
        double delta2 = data[i] - mean;
        M2 += delta*delta2;
    }
    std::vector<double> statistics = {mean, M2/count, M2/(count-1)};
    return statistics;
}

int main() {
    int N = 100;
    double s = 0.01;
    double ds = 0.0058;
    double R_sun = 6.955E8;
    
    std::random_device seeder;
    std::vector<std::mt19937_64> gens;
    for (int i = 0; i < omp_get_max_threads(); ++i) {
        std::mt19937_64 gen(seeder());
        gens.push_back(gen);
    }        
    
    std::vector<double> rs(N);
    std::vector<int> Steps(N);
    for (long i = 0; i < N; ++i) {
        photon p;
        p.scatter(gens);
        rs[i] = p.getR();
        Steps[i] = p.getSteps();
    }
    
    std::vector<double> r_avg = stats(rs);
    std::cout << "Average distance: " << r_avg[0] << " +/- " << std::sqrt(r_avg[2]) << std::endl;
    std::cout << "Number of steps: " << Steps[0] << ", " << Steps[3] << std::endl;
    return 0;
}
