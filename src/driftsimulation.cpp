#include <iostream>
#include <random>
#include <Rcpp.h>
#include "driftsimulation.h"

// constructor of the DriftSimulation class
// which is a discrete-time simulation of N haploid
// individuals
DriftSimulation::DriftSimulation(
        int const s_N // number of individuals
       ,double const s_v // value of winning
       ,double const s_c // cost of losing
       ,bool const s_pure // whether simulation works with pure or mixed strategy
       ,double const s_mu // mutation probability from H->D and D->H
       ,double const s_pHawk_init // initial frequency of hawks
       ,int const s_max_time // maximum duration of the simulation
       ,int s_output_nth_generation
       ,double const s_sd_pHawkMixed=0.0
       ) :
    rng_r((std::random_device())()) // init random number generator
    ,uniform(0.0, 1.0)
    ,N{s_N}
    ,v{s_v}
    ,c{s_c}
    ,is_pure{s_pure}
    ,payoff_matrix{{s_v/2.0, 0},{s_v, (s_v - s_c)/2.0}}
    ,mu{s_mu}
    ,max_time{s_max_time}
    ,pHawk_init{s_pHawk_init}
    ,output_nth_generation{s_output_nth_generation}
    ,sd_pHawkMixed{s_sd_pHawkMixed}
{
    if (output_nth_generation < 1)
    {
        std::ostringstream msg;

       msg << "Output is produced after " << output_nth_generation << " generations, which is invalid. A minimum value of 1 is required";

        throw msg.str();
    }
    else if (output_nth_generation >= max_time)
    {
        std::ostringstream msg;

        msg << "Output is produced after " << output_nth_generation << " generations, which is invalid (larger than maximum number of generations), setting to a minimum value of 1";

        throw msg.str();
    }

    double init_prob_hawk = is_pure ? 0.0 : pHawk_init;

    for (int pop_idx = 0; pop_idx < N; ++pop_idx)
    {
        Individual init_ind;

        init_ind.is_hawk = uniform(rng_r) < pHawk_init;
        init_ind.payoff = 0.0;

        // initialize the mixed probability of becoming a hawk
        init_ind.prob_hawk = init_prob_hawk;

        pop.push_back(init_ind);
    } // end for (int pop_idx)

} // end DriftSimulation::DriftSimulation()

// let hawks and doves interact with each other
// collect payoffs
// and then have them reproduce
void DriftSimulation::interact_reproduce()
{
    // shuffle the vector of individuals
    // so that we can make pairwise interactions
    // along the pop-vector's indices (i.e., 0 with 1,
    // 2 with 3, etc..) among random individuals
    std::shuffle(std::begin(pop),
                 std::end(pop),
                 rng_r
                 );

    bool ind1_is_hawk, ind2_is_hawk;
    // make baseline fitness variable to deal with neg numbers
    double baseline_fitness = c;

    std::vector <double> payoff_vector(pop.size());

    std::vector <Individual> next_generation(pop.size());

    // loop through individuals in pairs and calculate payoffs
    for (size_t ind_idx = 0; ind_idx < pop.size(); ind_idx+=2)
    {
        ind1_is_hawk = pop[ind_idx].is_hawk;
        ind2_is_hawk = pop[ind_idx + 1].is_hawk;

        // if we have HH we need to calculate
        // payoff stochastically
        if (ind1_is_hawk && ind2_is_hawk)
        {
            if (uniform(rng_r) < 0.5)
            {
                pop[ind_idx].payoff = baseline_fitness + v;
                pop[ind_idx + 1].payoff = baseline_fitness + -c;
            }
            else
            {
                pop[ind_idx].payoff = baseline_fitness + -c;
                pop[ind_idx + 1].payoff = baseline_fitness + v;
            }
        }
        else
        {
            // calculate payoffs for ind 1
            pop[ind_idx].payoff = baseline_fitness +
                payoff_matrix[ind1_is_hawk][ind2_is_hawk];

            // calculate payoffs for ind 2
            pop[ind_idx + 1].payoff = baseline_fitness +
                payoff_matrix[ind2_is_hawk][ind1_is_hawk];
        } // end if-else HH

        payoff_vector[ind_idx] = pop[ind_idx].payoff;

        assert(payoff_vector[ind_idx] >= 0.0);

        payoff_vector[ind_idx + 1] = pop[ind_idx + 1].payoff;

        assert(payoff_vector[ind_idx + 1] >= 0.0);

    } // end for int idx

    // make probability distribution based on payoffs
    std::discrete_distribution<int> parent_sampler(
        payoff_vector.begin()
        ,payoff_vector.end()
    );

    int sampled_parent;

    for (size_t ind_idx = 0; ind_idx < next_generation.size(); ++ind_idx)
    {
        sampled_parent = parent_sampler(rng_r);
        next_generation[ind_idx].is_hawk = pop[sampled_parent].is_hawk;
        next_generation[ind_idx].prob_hawk = pop[sampled_parent].prob_hawk;

        if (is_pure)
        {
            if (uniform(rng_r) < mu)
            {
                next_generation[ind_idx].is_hawk =
                    !next_generation[ind_idx].is_hawk;
            }
        }
        else
        {
            if (uniform(rng_r) < mu)
            {
                next_generation[ind_idx].prob_hawk =
                    next_generation[ind_idx].prob_hawk +
                        R::rnorm(0.0, sd_pHawkMixed);

                // set boundaries
                next_generation[ind_idx].prob_hawk =
                    next_generation[ind_idx].prob_hawk > 1.0 ?
                        1.0
                        :
                        (next_generation[ind_idx].prob_hawk < 0.0 ?
                            0.0
                            :
                            next_generation[ind_idx].prob_hawk);
            }

            // determine phenotype of this new individual
            next_generation[ind_idx].is_hawk = uniform(rng_r) <
                next_generation[ind_idx].prob_hawk;
        }
    }

    // replace current generation with the next
    pop = next_generation;

} // end DriftSimulation::interact_reproduce()


// get the stats
void DriftSimulation::write_data(
            int const generation_t
            ,Rcpp::NumericVector &freq_Hawk
            ,Rcpp::NumericVector &mean_pHawk
            ,Rcpp::NumericVector &sd_pHawk
            ,Rcpp::NumericVector & generation_vector)
{
    int n_hawk = 0;
    double p_hawk = 0.0;
    double ss_hawk = 0.0;
    double p;

    for (int ind_idx = 0; ind_idx < pop.size(); ++ind_idx)
    {
        n_hawk += pop[ind_idx].is_hawk;
        p = pop[ind_idx].prob_hawk;

        p_hawk += p;
        ss_hawk += p * p;
    }

    double mean_p_hawk = p_hawk/pop.size();

    freq_Hawk.push_back((double)n_hawk/pop.size());
    mean_pHawk.push_back(mean_p_hawk);
    sd_pHawk.push_back(ss_hawk/pop.size() - mean_p_hawk * mean_p_hawk);
    generation_vector.push_back(generation_t);
}

// run the actual simulation
Rcpp::DataFrame DriftSimulation::run()
{
    // calculate number of rows of the data frame
    int nrows = floor(static_cast<double>(max_time)/
                      output_nth_generation);

    // initialize the data
    Rcpp::NumericVector freq_Hawk{};
    Rcpp::NumericVector mean_pHawk{};
    Rcpp::NumericVector sd_pHawk{};
    Rcpp::NumericVector generation_vector{};



    for (int generation = 1;
         generation <= max_time;
         ++generation)
    {
        interact_reproduce();

        // write data every nth generation
        // to vector
        if (generation % output_nth_generation == 0)
        {
            write_data(generation
                            ,freq_Hawk
                           ,mean_pHawk
                           ,sd_pHawk
                           ,generation_vector);

            // check for user interruption
            Rcpp::checkUserInterrupt();
        }
    }

    // put everything in data.frame, in case
    // we want to collect different aspects of the dataset
    Rcpp::DataFrame simulation_result =
       Rcpp::DataFrame::create(
           Rcpp::Named("generation") = generation_vector
            ,Rcpp::Named("freq_Hawk") = freq_Hawk
            ,Rcpp::Named("mean_pHawkMixed") = mean_pHawk
            ,Rcpp::Named("sd_pHawkMixed") = sd_pHawk
        );

    return(simulation_result);
} // end DriftSimulation::run

