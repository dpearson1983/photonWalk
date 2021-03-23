#ifndef _PHOTON_HPP_
#define _PHOTON_HPP_

#include <vector>
#include <random>
#include <cmath>
#include <omp.h>

class photon{
    std::vector<double> pos;
    double r;
    int Steps;
        
public:
    photon() {
        this->pos = {0.0, 0.0, 0.0};
        this->Steps = 0;
    }
    
    void scatter(std::vector<std::mt19937_64> &gen) {
        std::vector<std::vector<double>> pos_i(omp_get_max_threads(), std::vector<double>(3));
        std::vector<int> N;
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            N.push_back(int(100000000/omp_get_max_threads() + 1));
            this->Steps += N[i];
        }
#pragma omp parallel
{
        int tid = omp_get_thread_num();
        std::uniform_real_distribution<double> uDist(0.0, 1.0);
        std::normal_distribution<double> nDist(0.01, 0.0058);
        for (int i = 0; i < N[tid]; ++i) {
            
            double s = nDist(gen[tid]);
            double theta = M_PI*uDist(gen[tid]);
            double phi = 2.0*M_PI*uDist(gen[tid]);
            double dx = s*std::cos(phi)*std::sin(theta);
            double dy = s*std::sin(phi)*std::sin(theta);
            double dz = s*std::cos(theta);
            pos_i[tid][0] += dx;
            pos_i[tid][1] += dy;
            pos_i[tid][2] += dz;
        }
}
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            this->pos[0] += pos_i[i][0];
            this->pos[1] += pos_i[i][1];
            this->pos[2] += pos_i[i][2];
        }
        this->r = std::sqrt(this->pos[0]*this->pos[0] + this->pos[1]*this->pos[1] + this->pos[2]*this->pos[2]);
    }
    
    std::vector<double> getPos() {
        return this->pos;
    }
    
    double getR() {
        return this->r;
    }
    
    int getSteps() {
        return this->Steps;
    }
};

#endif
