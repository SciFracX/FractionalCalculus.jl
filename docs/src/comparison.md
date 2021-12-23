# Comparison with other tools

Here we compare the current tools can be used to computing fractional calculus with FractionalCalculus.jl.

## Matlab

While most packages are encoded in Matlab, we would elaborate more on this section

### Matrix Discretization

The [Matrix Discretization method](https://www.mathworks.com/matlabcentral/fileexchange/22071-matrix-approach-to-discretization-of-odes-and-pdes-of-arbitrary-real-order) is a Matlab toolbox developed by Professor Igor Podlubny. FractionalCalculus.jl also support this method, the relating API is ```RLDiff_Matrix```.

### FOTF

[FOTF](https://www.mathworks.com/matlabcentral/fileexchange/60874-fotf-toolbox) is a Matlab toolbox developed by Professor Dingyu Xue, it can be used to modeling and analysis fractional order systems. FOTF contains various methods to compute fractional calculus, FractionalCalculus.jl adapted all the methods in this toolbox, the relating methods are:

```GL_High_Precision```

### Chebfun

It is noteworthy that [Chebfun](https://www.chebfun.org/) also can be used to compute fractional calculus.

[Detailed usage](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/23972/versions/22/previews/chebfun/examples/integro/html/FracCalc.html)

There are only one algorithm in Chebfun regarding Riemann Liouville fractional derivative and integral.

## Python



## R