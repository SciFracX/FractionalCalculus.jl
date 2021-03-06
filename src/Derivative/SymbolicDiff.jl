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

        @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys));
        # CONSTANTS AND POWERS
        @acrule sqrt(~x) => sqrt(pi)/2
        #@acrule((~x) => 2*sqrt(~x/pi))
        @acrule (~x)^(~α) => gamma(~α+1)/gamma(~α+1/2)*(~x)^(~α-1/2)

        # BINOMIALS
        @acrule 1/sqrt(~x) => 0
        @acrule sqrt(1+~x) => 1/sqrt(pi*~x) + atan(sqrt(~x))/sqrt(pi)
        @acrule sqrt(1-~x) => 1/sqrt(pi*~x) - atanh(sqrt(~x))/sqrt(pi)
        @acrule 1/sqrt(1+~x) => 1/(sqrt(pi*~x)*(1+~x))
        @acrule 1/sqrt(1-~x) => 1/(sqrt(pi*~x)*(1-~x))
        @acrule 1/(1+~x) => (sqrt(1+~x)-sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1+~x)^(3/2))
        @acrule 1/(1-~x) => (sqrt(1-~x)+sqrt(~x)*asinh(sqrt(~x)))/(sqrt(pi*~x)*(1-~x)^(3/2))
        @acrule 1/(1+~x)^(3/2) => (1-~x)/(sqrt(pi*~x)*(1+~x)^2)
        @acrule 1/(1-~x)^(3/2) => (1+~x)/(sqrt(pi*~x)*(1-~x)^2)
        @acrule 1/(1-~x)^(~p) => -Incomplete_beta(~x, -1/2, ~p+1/2)/(2*sqrt(pi)*(1-~x)^(~p+1/2))
        @acrule 1/sqrt(~x*(1-~x)) => (ellipe(~x)-(1-~x)*ellipk(~x))/(~x*(1-~x)*sqrt(pi))
        @acrule 1/sqrt(~x*(1+~x)) => (ellipe(~x/(1+~x))-ellipk(~x/(1+~x)))/(~x*sqrt(pi*(1+~x)))
        @acrule sqrt(~x/(1-~x)) => ellipe(~x)/((1-~x)*sqrt(pi))
        # Doesn't support Legendre
        @acrule 1/(sqrt(~x)*(1+~x)) => -sqrt(pi)/(2*(1+~x)^(3/2))
        @acrule 1/(sqrt(~x)*(1-~x)) => sqrt(pi)/(2*(1-~x)^(3/2))
        @acrule sqrt(~x)/~x => sqrt(pi)/(2*(1+~x)^(3/2))
        @acrule sqrt(~x)/(-~x) => sqrt(pi)/(2*(1-~x)^(3/2))
        @acrule sqrt(~x)/(1-~x)^2 => sqrt(pi)*(2+~x)/(4*(1-~x)^(5/2)) 
        #@acrule((~x)^(~p)/(1-(~x))^((~p)+3/2) => (((~p)+1/2+1/2*(~x))*gamma(~p+1)*(~x)^((~p)-1/2))/(gamma((~p+3/2)*(1-(~x))^(2+(~p))))
        @acrule sqrt((1-~x)/~x) => (ellipe(~x)-ellipk(~x))/(~x*sqrt(pi))
        @acrule sqrt((1+~x)/~x) => ((1+~x)*ellipe(~x/(1+~x)) - ellipk(~x/(1+~x)))/(~x*sqrt(pi*(1+~x)))
        @acrule  (~x)^(~p)/(1+~x)^~r => gamma(~p+1)/gamma(~p+1/2)*(~x)^(~p-1/2)* Hypergeometric1F1(~r, ~p+1, p+1/2, ~x)
        @acrule  (~x)^(~p)/(1-~x)^~r => gamma(~p+1)/gamma(~p+1/2)*(~x)^(~p-1/2)* Hypergeometric1F1(~r, ~p+1, p+1/2, -~x)



        # EXPONENTIAL AND RELATED FUNCTIONS
        @acrule exp(~x) => 1/sqrt(pi*~x) + exp(~x)*erf(sqrt(~x))
        @acrule exp(-~x) => 1/sqrt(pi*~x) - 2/sqrt(pi)*dawson(sqrt(-~x))
        @acrule exp(~x)*erf(sqrt(~x)) => exp(~x)
        @acrule dawson(sqrt(~x)) => 1/2*sqrt(pi)*exp(-~x)
        @acrule exp(~x)*erfc(sqrt(~x)) => 1/sqrt(pi*~x) - exp(~x)*erfc(sqrt(~x))
        @acrule exp(~x)*erfc(-sqrt(~x)) => 1/sqrt(pi*~x) + exp(~x)*erf(-sqrt(~x))
        @acrule erf(~x) => exp(-1/2*~x)*besseli(0, 1/2*~x)
        @acrule exp(~x)/sqrt(~x) => 1/2*sqrt(pi)*exp(1/2*~x)*(besseli(1, 1/2*~x)+besseli(0, 1/2*~x))
        @acrule exp(-~x)/sqrt(~x) => 1/2*sqrt(pi)*exp(-1/2*~x)*(besseli(1, 1/2*~x)-besseli(0, 1/2*~x)) 
        @acrule cosh(sqrt(~x))/sqrt(~x) => sqrt(pi/~x)*besseli(1, sqrt(~x))/2
        @acrule  (1-cos(sqrt(~x)))/~x => sqrt(pi)/2*Struve0(sqrt(~x)) 
        @acrule  (1-cosh(sqrt(~x)))/~x => sqrt(pi)/2*MStruve0(sqrt(~x)) 

        @acrule sin(~x) => sin(~x+pi/4)-sqrt(2)*(AuxiliaryFresnelCos(sqrt(2*~x/π)))
        @acrule cos(~x) => 1/(sqrt(π*~x))+cos(~x+1/4*π)-sqrt(2)*AuxiliaryFresnelSin(sqrt(2*~x/π))
        @acrule sinh(~x) => dawson(sqrt(~x))/sqrt(pi)-(exp(~x)*erf(sqrt(~x)))/2
        @acrule cosh(~x) => 1/sqrt(pi*~x)+exp(~x)*erf(sqrt(~x))/2-dawson(sqrt(~x))/sqrt(pi)
        @acrule asin(sqrt(~x))/sqrt(1-~x) => sqrt(pi)/(2*(1-~x))
        @acrule atan(sqrt(~x)) => 1/2*sqrt(pi/(1+~x))
        @acrule asinh(sqrt(~x))/sqrt(1+~x) => sqrt(pi)/(2*(1+~x))
        @acrule atanh(sqrt(~x)) => 1/2*sqrt(pi/(1-~x))

        # BESSEL AND STRUVE FUNCTIONS
        @acrule besselj(0, sqrt(~x)) => cos(sqrt(~x))/sqrt(pi*(~x))
        @acrule besselj(1, sqrt(~x))/sqrt(~x) => (cos(sqrt(~x))+sqrt(~x)*sin(sqrt(~x))-1)/(sqrt(~x)*(~x)^(3/2))
        @acrule  besselj(~v, sqrt(~x))/(~x)^(~v/2) => 1/(2^(~v)*gamma(~v+1)*sqrt(pi*~x))- Struve(~v+1/2, sqrt(~x))/(sqrt(2)*(~x)^(~v/2+1/4))
        @acrule sqrt(~x)*besselj(1, sqrt(~x)) => sin(sqrt(~x))/sqrt(pi)
        @acrule (~x)^(~v/2)*besselj(1, sqrt(~x)) => (~x)^((~v)/2-1/4)*besselj(~v-1/2, sqrt(~x))/sqrt(2)
        @acrule besseli(0, sqrt(~x)) => cosh(sqrt(~x))/sqrt(pi*~x)
        @acrule besseli(1, sqrt(~x))/sqrt(~x) => (cosh(sqrt(~x))-sqrt(~x)*sinh(sqrt(~x))-1)/(sqrt(pi)*(~x)^(3/2))
        @acrule  besseli(~v, sqrt(~x))/((~x)^(~v/2)) => 1/(2^~v*gamma(~v+1)*sqrt(pi*(~x))) + Legendre(~v+1/2, sqrt(~x))/(sqrt(2)*(~x)^(~v/2+1/4))
        @acrule sqrt(~x)*besseli(1, sqrt(~x)) => sinh(sqrt(~x))/sqrt(pi)
        @acrule  (~x)^(~v/2)*besseli(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*besseli(~v-1/2, sqrt(~x))/sqrt(2) 
        @acrule exp(-~x)*besseli(0, ~x) => exp(-2*~x)/sqrt(pi*(~x))
        @acrule  Struve(0, sqrt(~x)) => sin(sqrt(~x))/sqrt(pi*(~x))
        @acrule  Struve(0, sqrt(~x))/(~x) => (sin(sqrt(~x)) - SinIntegral(sqrt(~x)))/(sqrt(pi)*(~x)^(3/2)) 
        @acrule sqrt(~x)*Struve(1, sqrt(~x)) =>(1-cos(sqrt(~x)))/sqrt(pi)
        @acrule (~x)^(~v/2)*Struve(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*Struve(~v-1/2, sqrt(~x))/sqrt(2)

        @acrule MStruve(0, sqrt(~x)) => sinh(sqrt(~x))/sqrt(pi*(~x))
        @acrule MStruve(0, sqrt(~x))/(~x) => (sinh(sqrt(~x)) - HyperSinIntegral(sqrt(~x)))/(sqrt(pi)*(~x)^(3/2))
        @acrule sqrt(~x)*MStruve(1, sqrt(~x)) => (cosh(sqrt(~x)-1))/sqrt(pi)
        @acrule  (~x)^(~v/2)*MStruve(~v, sqrt(~x)) => (~x)^((~v)/2-1/4)*MStruve(~v-1/2, sqrt(~x))/sqrt(2)
        
        # MISCELLANEOUS FUNCTIONS
        @acrule log(~x) => log(4*~x)/sqrt(pi*~x)
        @acrule sqrt(~x)*log(~x) => sqrt(pi)/2*(log(~x/4)+2)
        @acrule log(~x)/sqrt(~x) => sqrt(pi)/~x

        @acrule ellipk(~x) => sqrt(pi)/(2*sqrt(~x*(1-~x)))
        @acrule ellipe(~x) => 1/2*sqrt(pi*(1-~x)/(~x))

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