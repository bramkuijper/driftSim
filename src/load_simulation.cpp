#include <Rcpp.h>
#include "driftsimulation.h"

//' Run a single replicate individual-based simulation
//' of a hawk-dove
//' game in a finite-sized population.
//' @param N Population size (integer).
//' @param v Value of winning a fight (floating-point value).
//' @param c Cost of losing a fight (floating-point value).
//' @param is_pure If \code{is_pure=TRUE}, Hawks and Doves are pure
//'    strategies, with mutation rates \code{mu} reflecting the probability
//'    that one pure strategy mutates into another. Alternatively, in case
//'    \code{is_pure=FALSE}, Hawks and Doves are the result of an evolving
//'    haploid locus that codes for probability \eqn{0 \leq pHawk \leq 1},
//'    which determines an individual's probability of playing Hawk.
//'    Here, \eqn{pHawk} gradually evolves. In case of a mutation event
//'    (ocurring at rate \code{mu}), the new value of the \eqn{pHwk}
//'    allele
//'    in individual \eqn{i} is given by \eqn{pHawk(i,t+1) = pHawk(i,t) +
//'    e}, where e is drawn from a normal distribution
//'    with mean 0 and standard deviation \code{sd_pHawkMixed} (see below)
//' @param mu Mutation rate from Hawk to Dove, from Dove to Hawk
//' @param max_time Maximum number of generations after which
//'    the simulation ends -- in case of extinction simulation ends
//'    earlier (integer).
//' @param pHawk_init Initial sampling probability of Hawks.
//'    If \code{pHawk_init = 0.5}, the initial proportion of Hawks should be
//'    close to one half. It is not exactly one half, as the initial number
//'    of Hawks is not deterministic, but is sampled from a binomial
//'    distribution
//' @param output_nth_generation Data is produced every
//'    \code{output_nth_generation} timesteps. Minimum number is 1 (output
//'    data every generation). This is done to save data in case you want to
//'    run the thing for
//'    say, values of \code{max_time = 50000}
//' @param sd_pHawkMixed In case the strategy is mixed, this is the standard
//'    deviation of the distribution from which new mutations in the probability
//'    of developing as a Hawk is sampled
//' @return A \code{data.frame()} that contains all kinds of statistics
//'    collected during the running of the simulation.
//' @examples
//' # quick test run of 10 generations and size of 10 individuals
//' # expecting ~100% hawks (but, yeah... sampling)
//' runSimulation(N=10, v=1.0, c=0.5, is_pure=TRUE, mu=0.01, max_time=10,
//'    pHawk_init=0.5, output_nth_generation=1, sd_pHawkMixed=0.0)
//'
//' # a much more extensive simulation of 1000 individuals for 10000 generations
//' # where the resulting data (every 10th generation) is stored in
//' # the variable sim.data
//' sim.data <- runSimulation(N=1000, v=1.0, c=2.0, is_pure=TRUE, mu=0.01,
//'    max_time=10000, pHawk_init=0.05, output_nth_generation=10,
//'    sd_pHawkMixed=0.0)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame runSimulation(
        int const N
        ,double const v
        ,double const c
        ,bool const is_pure
        ,double const mu
        ,int max_time
        ,double pHawk_init
        ,int output_nth_generation
        ,double const sd_pHawkMixed
        ) {

    DriftSimulation df =
        DriftSimulation(N
                ,v
                ,c
                ,is_pure
                ,mu
                ,pHawk_init
                ,max_time
                ,output_nth_generation
                ,sd_pHawkMixed);

    Rcpp::DataFrame output = df.run();

    return(output);
} //
