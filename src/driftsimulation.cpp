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
       ,double const s_mu // mutation probability from H->D and D->H
       ,double const s_pHawk_init // initial frequency of hawks
       ,int const s_max_time // maximum duration of the simulation
       ,int s_output_nth_generation) :
    rng_r((std::random_device())()) // init random number generator
    ,uniform(0.0, 1.0)
    ,N{s_N}
    ,v{s_v}
    ,c{s_c}
    ,payoff_matrix{{(s_v - s_c)/2.0, s_v},{0, s_v/2.0}}
    ,mu{s_mu}
    ,max_time{s_max_time}
    ,pHawk_init{s_pHawk_init}
    ,output_nth_generation{s_output_nth_generation}
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

    for (int pop_idx = 0; pop_idx < N; ++pop_idx)
    {
        Individual init_ind;

        init_ind.is_hawk = uniform(rng_r) < pHawk_init;
        init_ind.payoff = 0.0;

        pop.push_back(init_ind);
    } // end for (int pop_idx)

} // end DriftSimulation::DriftSimulation()

// let hawks and doves interact with each other and collect payoffs
int DriftSimulation::interact_reproduce()
{
    // shuffle the vector of individuals
    // so that we can make pairwise interactions
    // along the pop-vector's indices (i.e., 0 with 1,
    // 2 with 3, etc..)
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
    for (int ind_idx = 0; ind_idx < pop.size(); ind_idx+=2)
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
        payoff_vector[ind_idx + 1] = pop[ind_idx + 1].payoff;

    } // end for int idx

    // make probability distribution based on payoffs
    std::discrete_distribution<int> parent_sampler(
        payoff_vector.begin()
        ,payoff_vector.end()
    );

    int sampled_parent;

    // TODO
    for (int ind_idx = 0; ind_idx < next_gen.size(); ++ind_idx)
    {
        sampled_parent = parent_sampler(rng_r);
        next_gen[ind_idx].is_hawk = pop[sampled_parent].is_hawk;

        if (uniform(rng_r) < mu)
        {
            next_gen[ind_idx].is_hawk = !next_gen[ind_idx].is_hawk;
        }
    }

    // replace current generation with the next
    pop = next_gen;

} // end DriftSimulation::interact_reproduce()

Rcpp::DataFrame DriftSimulation::run()
{
    Rcpp::Rcout << "starting simulation..." << std::endl;

    // calculate number of rows of the data frame
    int nrows = floor(static_cast<double>(max_time)/
                      output_nth_generation);

    // initialize the data
    Rcpp::NumericVector pHawk(nrows);
    Rcpp::NumericVector unPaired(nrows);

    Rcpp::DataFrame simulation_result =
       Rcpp::DataFrame::create(
           Rcpp::Named("pHawk") = pHawk
            ,Rcpp::Named("unpaired") = unPaired
        );


    for (int generation = 1;
         generation <= max_time;
         ++generation)
         )
    {
        interact_reproduce();
    }

    return(simulation_result);
} // end DriftSimulation::run

