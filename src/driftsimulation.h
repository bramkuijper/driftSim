#ifndef DRIFT_SIMULATION_HPP
#define DRIFT_SIMULATION_HPP

#include <vector>
#include <random>
#include <Rcpp.h>

struct Individual {
    bool is_hawk;
    double payoff;
    double prob_hawk;
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
        bool is_pure;
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

        // standard deviation on the mixed strategy
        // evolution
        double sd_pHawkMixed;

        // let all hawks and doves interact
        void interact_reproduce();

        // append the data to the NumericVector data
        void write_data(
                int const generation_t
                ,Rcpp::NumericVector &freq_Hawk
                ,Rcpp::NumericVector &mean_pHawk
                ,Rcpp::NumericVector &sd_pHawk
                ,Rcpp::NumericVector & generation_vector);


    public:
       DriftSimulation(
               int const N
               ,double const v
               ,double const c
               ,bool const is_pure
               ,double const mu
               ,double const pHawk_init
               ,int const max_time
               ,int output_nth_generation
               ,double const sd_pHawkMixed
               );

       Rcpp::DataFrame run();
};

#endif
