#=
abstract type Symbolic end

abstract type AbstractVariable <: Symbolic end

const SymbolicVariable = Union{Symbol, AbstractVariable}
=#

incomplete_beta(z, a, b) = beta(a, b)*beta_inc(a, b, z)[1]

SEMIDIFFRULES = [
    # CONSTANTS AND POWERS
    @rule sqrt(~x) => sqrt(pi)/2
    @rule (~x) => 2*sqrt(~x/pi)
    @rule (~x)^(~α) => gamma(~α+1)/gamma(~α+1/2)*(~x)^(~α-1/2)

    # BINOMIALS
    @rule sqrt(1+~x) => 1/sqrt(pi*~x) + atan(sqrt(~x))/sqrt(pi)
    @rule sqrt(1-~x) => 1/sqrt(pi*~x) - atanh(sqrt(~x))/sqrt(pi)
    @rule 1/sqrt(1+~x) => 1/(sqrt(pi*~x)*(1+~x))
    @rule 1/sqrt(1-~x) => 1/(sqrt(pi*~x)*(1-~x))
    @rule 1/(1+~x) => (sqrt(1+~x)-sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1+~x)^(3/2))
    @rule 1/(1-~x) => (sqrt(1-~x)+sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1-~x)^(3/2))
    @rule 1/(1+~x) => (1-~x)/(sqrt(pi*~x)*(1+~x)^2)
    @rule 1/(1-~x) => (1+~x)/(sqrt(pi*~x)*(1-~x)^2)
    @rule 1/(1-~x)^(~p) => -incomplete_beta(~x, -1/2, ~p+1/2)/(2*sqrt(pi)*(1-~x)^(~p+1/2))
    #FIXME: Elliptic integral functions...

    # EXPONENTIAL AND RELATED FUNCTIONS
    @rule exp(~x) => 1/sqrt(pi*~x) + exp(~x)*erf(sqrt(~x))
    @rule exp(-~x) => 1/sqrt(pi*~x) - 2/sqrt(pi)*dawson(sqrt(-~x))
    @rule exp(~x)*erf(sqrt(~x)) => exp(~x)
    @rule dawson(sqrt(~x)) => 1/2*sqrt(pi)*exp(-~x)
    @rule exp(~x)*erf(sqrt(~x)) => 1/sqrt(pi*~x) - exp(~x)*erf(sqrt(~x))
    @rule exp(~x)*erf(-sqrt(~x)) => 1/sqrt(pi*~x) + exp(~x)*erf(-sqrt(~x))
    @rule erf(~x) => exp(-1/2*~x)*besseli(0, 1/2*~x)
]