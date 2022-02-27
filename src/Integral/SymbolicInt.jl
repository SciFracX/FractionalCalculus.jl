
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
        @acrule(exp(~x)*erf(sqrt(~x)) => exp(~x)-1)
        @acrule(dawson(~x) => 1/2*sqrt(pi)*(1-exp(-~x)))
        @acrule(exp(~x)*erfc(sqrt(~x)) => 1-exp(~x)*erfc(sqrt(~x)))
        @acrule(exp(~x)*erfc(-sqrt(~x)) => exp(~x)*erfc(sqrt(~x))-1)
        @acrule(erf(sqrt(~x)) => ~x*exp(-1/2*~x)*(besseli(1, 1/2*~x)+besseli(0, 1/2*~x)))
        @acrule(exp(~x)/sqrt(~x) => sqrt(pi)*exp(1/2*~x)*besseli(0, 1/2*~x))
        @acrule(exp(-~x)/sqrt(~x) => sqrt(pi)*exp(-1/2*~x)*besseli(0, 1/2*~x))
        @acrule(cosh(sqrt(~x)) => sqrt(pi)*besseli(0, sqrt(~x)))
        # Struve function and Modified Struve function @acrule (1-cos(sqrt(~x)))/~x => sqrt(pi)/2*
        @acrule(sinh(~x) => (exp(~x)*erf(sqrt(~x)))/2-dawson(sqrt(~x))/sqrt(pi))
        @acrule(cosh(~x) => (exp(~x)*erf(sqrt(~x)))/2+dawson(sqrt(~x))/sqrt(pi))
        @acrule(asin(sqrt(~x))/sqrt(1-~x) => -sqrt(pi)/2*log(1-~x))
        @acrule(atan(sqrt(~x)) => sqrt(pi)*(sqrt(1+~x)-1))
        @acrule(asinh(sqrt(~x))/sqrt(1+~x) => sqrt(pi)/2*log(1+~x))
        @acrule(atanh(sqrt(~x)) => sqrt(pi)*(1-sqrt(1-~x)))

        # BESSEL AND STRUVE FUNCTIONS
        @acrule(besselj(0, sqrt(~x)) => sin(sqrt(~x))/sqrt(pi*(~x)))
        @acrule(besselj(1, sqrt(~x)) => 2*(1-cos(sqrt(~x)))/sqrt(pi*~x))
        # Struve function doesn't support
        @acrule(sqrt(~x)*besselj(1, sqrt(~x)) => (2*sin(sqrt(~x))-2*sqrt(~x)*cos(sqrt(~x)))/sqrt(pi))
        @acrule((~x)^(~v/2)*besselj(sqrt(~x)) => (~x)^((~v)/2+1/4)*besselj(~v+1/2, sqrt(~x)))
        @acrule(besseli(0, sqrt(~x)) => sinh(sqrt(~x))/sqrt(pi*~x))
        @acrule(besseli(1, sqrt(~x))/sqrt(~x) => 2*(1-cosh(sqrt(~x)))/sqrt(pi*~x))
        @acrule(sqrt(~x)*besseli(1, sqrt(~x)) => (2*sqrt(~x)*cosh(sqrt(~x))-2*sinh(sqrt(~x)))/sqrt(pi))
        @acrule(exp(-~x)*besseli(0, ~x) => erf(sqrt(2*~x))/sqrt(2))

        # MISCELLANEOUS FUNCTIONS
        @acrule(log(~x) => 2*sqrt(~x/pi)*(log(4*~x-2)))
        @acrule(sqrt(~x)*log(~x) => (~x*sqrt(pi))/2*(log(~x/4)+1))
        @acrule(log(~x)/sqrt(~x) => sqrt(pi)*log(~x/4))
]
#=
EXP_RELATED_FUN_RULES = [

]

BESSEL_STRUVE_FUN_RULES = [

]

MISCELLANEOUS_FUN_RULES = [

]
=#

INTRULES = Chain(SEMIINTRULES)

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
semiint(x) = INTRULES(x)