#include <Rcpp.h>
#include "driftsimulation.h"

// ' Run a single replicate individual-based simulation
// ' of a hawk-dove
// ' game in a finite-sized population.
// ' @param N Population size (integer).
// ' @param v A Value of winning a fight (floating-point value).
// ' @param c A Cost of losing a fight (floating-point value).
// ' @param max_time Maximum number of generations after which
// '    the simulation ends in case there is no extinction (integer).
// ' @param pHawk_init Initial sampling probability of Hawks. It depends
// '    on the realization of the random process how many Hawks actually
// '    will be there (integer).
// ' @param output_nth_generation Do you want to obtain data of every
// '    generation or will there be an interval of
// '    \code{output_nth_generation} timesteps. Minimum number is 1.
// ' @return A \code{data.frame()} that contains all kinds of statistics
// '    collected during the running of the simulation.
// ' @examples
// '
// ' @export
// [[Rcpp::export]]
Rcpp::DataFrame runSimulation(
        int const N
        ,double const v
        ,double const c
        ,double const mu
        ,int max_time
        ,double pHawk_init
        ,int output_nth_generation
        ) {

    DriftSimulation df =
        DriftSimulation(N
                ,v
                ,c
                ,mu
                ,pHawk_init
                ,max_time
                ,output_nth_generation);

    Rcpp::DataFrame output = df.run();

    return(output);
} //
