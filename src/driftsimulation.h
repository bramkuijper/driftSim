#ifndef DRIFT_SIMULATION_HPP
#define DRIFT_SIMULATION_HPP

#include <vector>
#include <random>
#include <Rcpp.h>

struct Individual {
    bool is_hawk;
    double payoff;
};

class DriftSimulation
{
    private:
        std::mt19937 rng_r; // random number generator
        std::uniform_real_distribution<double> uniform;
        std::vector <Individual> pop;

        int N;
        double v;
        double c;
        double payoff_matrix[2][2];
        double mu;
        int max_time;
        double pHawk_init;

        // output produced every nth generation
        // if 1, output is produced every generation
        // if 10, output is produced every 10 generations
        // etc
        //
        // 0 will be ignored and set at 1
        int output_nth_generation;

        // let all hawks and doves interact
        void interact_reproduce();

        // append the data to the NumericVector data
        void write_data(Rcpp::NumericVector &data);


    public:
       DriftSimulation(
               int const N
               ,double const v
               ,double const c
               ,double const mu
               ,double const pHawk_init
               ,int const max_time
               ,int output_nth_generation
               );

       Rcpp::DataFrame run();
};

#endif
