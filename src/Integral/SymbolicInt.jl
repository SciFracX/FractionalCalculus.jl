
SEMIINTRULES = [
        # CONSTANTS AND POWERS
        @acrule(sqrt(~x) => sqrt(pi)/2*~x)
        #@acrule((~x) => 4*~x^(3/2)/(3*sqrt(pi)))
        @acrule((~x)^(~α) => gamma(~α+1)/gamma(~α+3/2)*(~x)^(~α+1/2))

        # BINOMIALS
        @acrule(sqrt(1+~x) => sqrt(~x/pi) + ((1+~x)*atan(sqrt(~x)))/sqrt(pi))
        @acrule(sqrt(1-~x) => sqrt(~x/pi) + ((1-~x)*atanh(sqrt(~x)))/sqrt(pi))
        @acrule(1/sqrt(1+~x) => 2/sqrt(pi)*atan(sqrt(~x)))
        @acrule(1/sqrt(1-~x) => 2/sqrt(pi)*atanh(sqrt(~x)))
        @acrule(1/(1+~x) => 2*asinh(sqrt(~x))/sqrt(pi*(1+~x)))
        @acrule(1/(1-~x) => 2*asin(sqrt(~x))/sqrt(pi*(1-~x)))
        @acrule(1/(1+~x)^(3/2) => 2/(1+~x)*sqrt(~x/pi))
        @acrule(1/(1-~x)^(3/2) => 2/(1-~x)*sqrt(~x/pi))
        @acrule(1/(1-~x)^(~p) => -incomplete_beta(~x, 1/2, ~p-1/2)/(sqrt(pi)*(1-~x)^(~p-1/2)))
        @acrule(1/sqrt(~x*(1-~x)) => 2*ellipk(~x)/sqrt(pi))
        @acrule(1/sqrt(~x*(1+~x)) => 2*ellipk(~x/(1+~x))/sqrt(pi*(1+~x)))
        @acrule(sqrt(~x/(1-~x)) => 2/sqrt(pi)*(ellipk(~x)-ellipe(~x)))
        @acrule(1/(sqrt(~x)*(1+~x)) => sqrt(pi/(1+~x)))
        @acrule(1/(sqrt(~x)*(1-~x)) => sqrt(pi/(1-~x)))
        @acrule(sqrt(~x)/~x => -sqrt(pi/(1+~x))+sqrt(pi))
        @acrule(sqrt(~x)/(-~x) => sqrt(pi/(1-~x))-sqrt(pi))
        @acrule(sqrt(~x)/(1-~x)^2 => sqrt(pi)*~x/(2*(1-~x)^(3/2)))
        @acrule(sqrt((1-~x)/~x) => 2/sqrt(pi)*ellipe(~x))
        @acrule(sqrt((1+~x)/~x) => 2*sqrt((1+~x)/pi)*ellipe(~x/(1+~x)))
        #Gauss hypergeometric functions

        # EXPONENTIAL AND RELATED FUNCTIONS
        @acrule(exp(~x) => exp(~x)*erf(sqrt(~x)))
        @acrule(exp(-~x) => 2/sqrt(pi)*dawson(sqrt(~x)))
]
#=
EXP_RELATED_FUN_RULES = [
        # EXPONENTIAL AND RELATED FUNCTIONS
]

BESSEL_STRUVE_FUN_RULES = [

]

MISCELLANEOUS_FUN_RULES = [

]
=#

RULES = Chain(SEMIINTRULES)

"""
# Symbolic fractional integral

    semiint(fun)

```semiint``` uses SymbolicUtils.jl to compute symbolic fractional integral.

### Example

```julia-repl
julia> semiint(x^4)
0.45851597901024005(x^4.5)
```
"""
semiint(x) = RULES(x)