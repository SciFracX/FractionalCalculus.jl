using SymbolicUtils, SpecialFunctions
using SymbolicUtils.Rewriters
using Symbolics

Incomplete_beta(z, a, b) = beta(a, b)*beta_inc(a, b, z)[1]

# These are direct computing, we need to change to the numerical methods in http://www.mymathlib.com/c_source/functions/fresnel_sin_cos_integrals/fresnel_auxiliary_sine_integral.c
AuxiliaryFresnelSin(x) = sqrt(2/pi)*quadgk(t->exp(-2x*t)*sin(t^2), 0, Inf)[1]
AuxiliaryFresnelCos(x) = sqrt(2/pi)*quadgk(t->exp(-2x*t)*cos(t^2), 0, Inf)[1]

@register_symbolic Incomplete_beta(z::Number, a::Number, b::Number)
@register_symbolic AuxiliaryFresnelCos(x::Number)
@register_symbolic AuxiliaryFresnelSin(x::Number)
@register_symbolic Hypergeometric1F1(a::Number, b::Number, c::Number, z::Number)
@register_symbolic Struve(p::Number, x::Number)# Here p is the order of Struve function
@register_symbolic MStruve(p::Number, x::Number)# Here p is the order of modified Struve function
@register_symbolic Legendre(p::Number, x::Number)# Here p is the order of Legendre function
@register_symbolic SinIntegral(x::Number)
@register_symbolic HyperSinIntegral(x::Number)



SEMIDIFFRULES = [

        @rule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys));
        # CONSTANTS AND POWERS
        #@rule sqrt(~x) => sqrt(pi)/2
        #@rule((~x) => 2*sqrt(~x/pi))
        @rule (~x)^(~α) => gamma(~α+1)/gamma(~α+1/2)*(~x)^(~α-1/2)

        # BINOMIALS
        #@rule 1/sqrt(~x) => 0
        @rule sqrt(1+~x) => 1/sqrt(pi*~x) + atan(sqrt(~x))/sqrt(pi)
        @rule sqrt(1-~x) => 1/sqrt(pi*~x) - atanh(sqrt(~x))/sqrt(pi)
        @rule 1/sqrt(1+~x) => 1/(sqrt(pi*~x)*(1+~x))
        @rule 1/sqrt(1-~x) => 1/(sqrt(pi*~x)*(1-~x))
        @rule 1/(1+~x) => (sqrt(1+~x)-sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1+~x)^(3/2))
        @rule 1/(1-~x) => (sqrt(1-~x)+sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1-~x)^(3/2))
        @rule 1/(1+~x)^(3/2) => (1-~x)/(sqrt(pi*~x)*(1+~x)^2)
        @rule 1/(1-~x)^(3/2) => (1+~x)/(sqrt(pi*~x)*(1-~x)^2)
        @rule 1/(1-~x)^(~p) => -Incomplete_beta(~x, -1/2, ~p+1/2)/(2*sqrt(pi)*(1-~x)^(~p+1/2))
        @rule 1/sqrt(~x*(1-~x)) => (ellipe(~x)-(1-~x)*ellipk(~x))/(~x*(1-~x)*sqrt(pi))
        @rule 1/sqrt(~x*(1+~x)) => (ellipe(~x/(1+~x))-ellipk(~x/(1+~x)))/(~x*sqrt(pi*(1+~x)))
        @rule sqrt(~x/(1-~x)) => ellipe(~x)/((1-~x)*sqrt(pi))
        # Doesn't support Legendre
        @rule 1/(sqrt(~x)*(1+~x)) => -sqrt(pi)/(2*(1+~x)^(3/2))
        @rule 1/(sqrt(~x)*(1-~x)) => sqrt(pi)/(2*(1-~x)^(3/2))
        @rule sqrt(~x)/~x => sqrt(pi)/(2*(1+~x)^(3/2))
        @rule sqrt(~x)/(-~x) => sqrt(pi)/(2*(1-~x)^(3/2))
        @rule sqrt(~x)/(1-~x)^2 => sqrt(pi)*(2+~x)/(4*(1-~x)^(5/2)) 
        #@rule((~x)^(~p)/(1-(~x))^((~p)+3/2) => (((~p)+1/2+1/2*(~x))*gamma(~p+1)*(~x)^((~p)-1/2))/(gamma((~p+3/2)*(1-(~x))^(2+(~p))))
        @rule sqrt((1-~x)/~x) => (ellipe(~x)-ellipk(~x))/(~x*sqrt(pi))
        @rule sqrt((1+~x)/~x) => ((1+~x)*ellipe(~x/(1+~x)) - ellipk(~x/(1+~x)))/(~x*sqrt(pi*(1+~x)))
        @rule  (~x)^(~p)/(1+~x)^~r => gamma(~p+1)/gamma(~p+1/2)*(~x)^(~p-1/2)* Hypergeometric1F1(~r, ~p+1, p+1/2, ~x)
        @rule  (~x)^(~p)/(1-~x)^~r => gamma(~p+1)/gamma(~p+1/2)*(~x)^(~p-1/2)* Hypergeometric1F1(~r, ~p+1, p+1/2, -~x)



        # EXPONENTIAL AND RELATED FUNCTIONS
        @rule exp(~x) => 1/sqrt(pi*~x) + exp(~x)*erf(sqrt(~x))
        @rule exp(-~x) => 1/sqrt(pi*~x) - 2/sqrt(pi)*dawson(sqrt(-~x))
        @rule exp(~x)*erf(sqrt(~x)) => exp(~x)
        @rule dawson(sqrt(~x)) => 1/2*sqrt(pi)*exp(-~x)
        @rule exp(~x)*erfc(sqrt(~x)) => 1/sqrt(pi*~x) - exp(~x)*erfc(sqrt(~x))
        @rule exp(~x)*erfc(-sqrt(~x)) => 1/sqrt(pi*~x) + exp(~x)*erf(-sqrt(~x))
        @rule erf(~x) => exp(-1/2*~x)*besseli(0, 1/2*~x)
        @rule exp(~x)/sqrt(~x) => 1/2*sqrt(pi)*exp(1/2*~x)*(besseli(1, 1/2*~x)+besseli(0, 1/2*~x))
        @rule exp(-~x)/sqrt(~x) => 1/2*sqrt(pi)*exp(-1/2*~x)*(besseli(1, 1/2*~x)-besseli(0, 1/2*~x)) 
        @rule cosh(sqrt(~x))/sqrt(~x) => sqrt(pi/~x)*besseli(1, sqrt(~x))/2
        @rule  (1-cos(sqrt(~x)))/~x => sqrt(pi)/2*Struve0(sqrt(~x)) 
        @rule  (1-cosh(sqrt(~x)))/~x => sqrt(pi)/2*MStruve0(sqrt(~x)) 

        @rule sin(~x) => sin(~x+pi/4)-sqrt(2)*(AuxiliaryFresnelCos(sqrt(2*~x/π)))
        @rule cos(~x) => 1/(sqrt(π*~x))+cos(~x+1/4*π)-sqrt(2)*AuxiliaryFresnelSin(sqrt(2*~x/π))
        @rule sinh(~x) => dawson(sqrt(~x))/sqrt(pi)-(exp(~x)*erf(sqrt(~x)))/2
        @rule cosh(~x) => 1/sqrt(pi*~x)+exp(~x)*erf(sqrt(~x))/2-dawson(sqrt(~x))/sqrt(pi)
        @rule asin(sqrt(~x))/sqrt(1-~x) => sqrt(pi)/(2*(1-~x))
        @rule atan(sqrt(~x)) => 1/2*sqrt(pi/(1+~x))
        @rule asinh(sqrt(~x))/sqrt(1+~x) => sqrt(pi)/(2*(1+~x))
        @rule atanh(sqrt(~x)) => 1/2*sqrt(pi/(1-~x))

        # BESSEL AND STRUVE FUNCTIONS
        @rule besselj(0, sqrt(~x)) => cos(sqrt(~x))/sqrt(pi*(~x))
        @rule besselj(1, sqrt(~x))/sqrt(~x) => (cos(sqrt(~x))+sqrt(~x)*sin(sqrt(~x))-1)/(sqrt(~x)*(~x)^(3/2))
        @rule  besselj(~v, sqrt(~x))/(~x)^(~v/2) => 1/(2^(~v)*gamma(~v+1)*sqrt(pi*~x))- Struve(~v+1/2, sqrt(~x))/(sqrt(2)*(~x)^(~v/2+1/4))
        @rule sqrt(~x)*besselj(1, sqrt(~x)) => sin(sqrt(~x))/sqrt(pi)
        @rule (~x)^(~v/2)*besselj(1, sqrt(~x)) => (~x)^((~v)/2-1/4)*besselj(~v-1/2, sqrt(~x))/sqrt(2)
        @rule besseli(0, sqrt(~x)) => cosh(sqrt(~x))/sqrt(pi*~x)
        @rule besseli(1, sqrt(~x))/sqrt(~x) => (cosh(sqrt(~x))-sqrt(~x)*sinh(sqrt(~x))-1)/(sqrt(pi)*(~x)^(3/2))
        @rule  besseli(~v, sqrt(~x))/((~x)^(~v/2)) => 1/(2^~v*gamma(~v+1)*sqrt(pi*(~x))) + Legendre(~v+1/2, sqrt(~x))/(sqrt(2)*(~x)^(~v/2+1/4))
        @rule sqrt(~x)*besseli(1, sqrt(~x)) => sinh(sqrt(~x))/sqrt(pi)
        @rule  (~x)^(~v/2)*besseli(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*besseli(~v-1/2, sqrt(~x))/sqrt(2) 
        @rule exp(-~x)*besseli(0, ~x) => exp(-2*~x)/sqrt(pi*(~x))
        @rule  Struve(0, sqrt(~x)) => sin(sqrt(~x))/sqrt(pi*(~x))
        @rule  Struve(0, sqrt(~x))/(~x) => (sin(sqrt(~x)) - SinIntegral(sqrt(~x)))/(sqrt(pi)*(~x)^(3/2)) 
        @rule sqrt(~x)*Struve(1, sqrt(~x)) =>(1-cos(sqrt(~x)))/sqrt(pi)
        @rule (~x)^(~v/2)*Struve(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*Struve(~v-1/2, sqrt(~x))/sqrt(2)

        @rule MStruve(0, sqrt(~x)) => sinh(sqrt(~x))/sqrt(pi*(~x))
        @rule MStruve(0, sqrt(~x))/(~x) => (sinh(sqrt(~x)) - HyperSinIntegral(sqrt(~x)))/(sqrt(pi)*(~x)^(3/2))
        @rule sqrt(~x)*MStruve(1, sqrt(~x)) => (cosh(sqrt(~x)-1))/sqrt(pi)
        @rule  (~x)^(~v/2)*MStruve(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*MStruve(~v-1/2, sqrt(~x))/sqrt(2)
        
        # MISCELLANEOUS FUNCTIONS
        @rule log(~x) => log(4*~x)/sqrt(pi*~x)
        @rule sqrt(~x)*log(~x) => sqrt(pi)/2*(log(~x/4)+2)
        @rule log(~x)/sqrt(~x) => sqrt(pi)/(~x)

        @rule ellipk(~x) => sqrt(pi)/(2*sqrt(~x*(1-~x)))
        @rule ellipe(~x) => 1/2*sqrt(pi*(1-~x)/(~x))

]

DIFFRULES = Chain(SEMIDIFFRULES)

"""
# Symbolic fractional differentiation

    semidiff(fun)

```semidiff``` uses SymbolicUtils.jl and Symbolics.jl to compute symbolic fractional differentiation.

### Example

```julia-repl
julia> using SymbolicUtils
julia> @syms x
julia> semidiff(log(x))
log(4x) / sqrt(πx)
```
"""
semidiff(x) = DIFFRULES(x)