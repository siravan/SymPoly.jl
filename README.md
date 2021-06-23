# SymPoly.jl

**SymPoly.jl** is a symbolic polynomial manipulation package with two main goals:

1. **SymPoly.jl** is a bridge between the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions and [Julia Algebra](https://github.com/JuliaAlgebra) polynomials (*DynamicPolynomials* in https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

2. **SymPoly.jl** provides a collections of advances *univariate* polynomial algorithms, including roof finding, extended greatest common divisor (GCD), polynomial factorization, and partial rational decomposition.

## General Polynomials

### Conversion between Symbolics Expressions and DynamicPolynomials

We can convert a Symbolics.jl expression to a DynamicPolynomials polynomial using `poly(eq, x => y)` function, where `eq` is a Symbolics expression, `x` is the name of the Symbolics variable and `y` is the name of `Polyvar` in the DynamicPolynomials. The conversion back from a DynamicPolynomials to a Symbolic expression is done with the help of `sym(p::AbstractPolynomial, y => x)` function.

```julia
julia> using SymPoly

julia> using SymbolicUtils

julia> using DynamicPolynomials

julia> @syms 𝑥
(𝑥,)

julia> @polyvar x
(x,)

julia> eq = 𝑥^2 - 4𝑥 + 3
3 + 𝑥^2 - (4𝑥)

julia> typeof(eq)
SymbolicUtils.Add{Number, Int64, Dict{Any, Number}, Nothing}

julia> p = poly(eq, 𝑥 => x)
x² - 4x + 3

julia> typeof(p)
Polynomial{true, Int64}

julia> sym(p, x => 𝑥)
3 + 𝑥^2 - (4𝑥)
```

Note that both Symbolics and DynamicPolynomials support multivariate polynomials; however, the transformation between the two is currently limited to univariate polynomials.

`poly` and `sym` have other forms. If the Symbolics expression has only one variable, we can simplify the second argument to just the name of the target variable (`poly(eq, 𝑥)` above). We can even skip the second argument altogether. It is case, a default name of `𝑦` is assigned to the second variable.

The following helper functions can be used to get various components of a polynomial expression or convert it to different standard forms:

* `var(p::AbstractPolynomial)`: returns the main polynomial variable.

* `deg(p::AbstractPolynomial)`: returns the degree of a polynomial (in the main variable).

* `leading(p::AbstractPolynomial)`: returns the leading coefficient (the coefficient of the highest non-zero power of `x`).

* `cont(p::AbstractPolynomial)`: returns the *content* of a polynomial with integer or rational coefficients, i.e., the GCD of the coefficients.

* `prim(p::AbstractPolynomial)`: returns a *primitive* polynomial, i.e., `p / cont(p)`, so the coefficients are co-primes.

### Basic Operations

The output of `poly` is a standard DynamicPolynomials and therefore supports all the basic operations on polynomials, including `+`, `-`, `*`, `÷`, `%`, and `^`. In addition, polynomial evaluation for a given value of the main variable, say `x₀`, is defined as `p(x₀)`.

### Greatest Common Divisor (GCD)

Finding the GCD of two polynomials is one of the basic blocks of many advanced polynomial manipulations (such as factorization). [Julia Algebra](https://github.com/JuliaAlgebra) provides a `gcd` function. **SymPoly** adds `gcdx`, the extended Euclid algorithm, which returns three values `g, s, t`, such as `g` is the GCD and `g = s*p + t*q`. Lets' look at an example:

```julia
julia> p = (x-1) * (x-4)
x² - 5x + 4

julia> q = (x-1) * (x+2)
x² + x - 2

julia> gcd(p, q)
-x + 1

julia> g, s, t = gcdx(p, q)
(x - 1//1, -1//6, 1//6)

julia> isequal(s*p + t*q, g)
true

```

`gcd(p, q)` is `-x + 1` instead of the expected `x - 1`, as `gcd` has a multiplicative ambiguity. Note that `gcdx` uses rational coefficients (in fact, `Rational{BigInt}`). This is a feature of most of the advanced algorithms in **SymPoly**.

### Root-Finding

We can find the real and complex roots of a polynomial with constant coefficients using `find_roots(p::AbstractPolynomial, x)`. It returns two arrays, the first one is a list of the real roots and the second one is a list of the complex ones. For example,

```julia
julia> p = (x-1) * (x-2) * (x^2+2x+5)
x⁴ - x³ + x² - 11x + 10

julia> r, s = find_roots(p, x)
([1.0, 2.0], Complex[-1.0 + 2.0im, -1.0 - 2.0im])
```

### Polynomial Factorization

[Polynomial factorization](https://en.wikipedia.org/wiki/Factorization_of_polynomials) is a very important operation and is required for symbolic integration, among other applications. From this point on, we are limited to polynomials with constant coefficients. Specially, we work with rational polynomials, whose coefficients are of type `Rational{BigInt}`. We can use the helper function `rationalize(p::AbstractPolynomial)` to convert a polynomial to a rational one. Currently, three factorization algorithms are provided:

1. `factor_schubert_kronecker(p::AbstractPolynomial)`, which as its name implies uses the Schubert-Kronecker algorithm. Note that this algorithm is based on old ideas (Friedrich Theodor von Schubert died in 1825 and Leopold Kronecker died in 1891) and is not efficient. However, it is usable for low-degree polynomials.

2. `decompose(p::AbstractPolynomial)` that uses the [Yun's algorithm](https://en.wikipedia.org/wiki/Square-free_polynomial) to decompose a polynomial into a list of co-prime and square-free factors. This is a fast algorithm but the decomposition is incomplete. For some applications, such as simplifying rational expressions, this is all that is needed. In other occasions, it can serve as the first step in a completet factorization (see below).

3. `factor(p::AbstractPolynomial)` is the main factorization entry point and currently uses a combination of square-free decomposition (first step) and the Schubert-Kronecker algorithm (second-step) to factor a polynomial. The plan is to switch to an efficient *Cantor–Zassenhaus algorithm* for the second step.

All three functions returns a value of type `FactoredPoly`, which is a list of `factor => power` pairs. We can access the factors and powers with the help of `factors`, `power`, `poly`, and the overloaded index operator (see below). Let's look at an example.

```julia
julia> p = (x-2) * (x+5)^2 * (x^2 + 2x + 5)
x⁵ + 10x⁴ + 26x³ - 75x - 250

julia> f = factor(p)
FactoredPoly(Pair{Any, Int64}[x - 2 => 1, x² + 2x + 5 => 1, x + 5 => 2])

julia> factors(f)
3-element Vector{Pair{Any, Int64}}:
       x - 2 => 1
 x² + 2x + 5 => 1
       x + 5 => 2

julia> f[3]
x + 5

julia> power(f, 3)
2

julia> isequal(f[1]^power(f,1) * f[2]^power(f,2) * f[3]^power(f,3), p)
true

julia> poly(f)
x⁵ + 10x⁴ + 26x³ - 75x - 250
```

### Partial Fraction Decomposition

[Partial Fraction Decomposition](https://en.wikipedia.org/wiki/Partial_fraction_decomposition) is the conversion of a fraction, where both the numerator and denominator are polynomials, to a sum of simpler fractions. Partial fraction decomposition has many applications, including in symbolic integration.

**SymPoly** uses the Hermite's method based on the square-free decomposition of the denominator to perform partial fraction decomposition. Let `p` be the numerator and `q` the denominator, where both are polynomials on the same variable. We can calculate the partial fraction decomposition of `p / q` by using `expand_frac(p, q)`. For example,

```julia
julia> p = x^4 + 6x^3 + 7x^2 +6x + 4
x⁴ + 6x³ + 7x² + 6x + 4

julia> q = x^6 + 2x^5 + 4x^4 + 4x^3 + 4x^2 + 2x + 1
x⁶ + 2x⁵ + 4x⁴ + 4x³ + 4x² + 2x + 1

julia> f = factor(p / q)
(2) / (x^2 + 1) + (-x^2 + 2*x + 2) / (x^2 + x + 1)^2

julia> factors(f)
2-element Vector{Pair{Any, Int64}}:
                (2) / (x² + 1) => 1
 (-x² + 2x + 2) / (x² + x + 1) => 2

julia> f[1]
(2) / (x² + 1)

julia> f[2]
(-x² + 2x + 2) / (x² + x + 1)

julia> a = numerator(f[1]) / denominator(f[1])^power(f,1)
(2) / (x² + 1)

julia> b = numerator(f[2]) / denominator(f[2])^power(f,2)
(-x² + 2x + 2) / (x⁴ + 2x³ + 3x² + 2x + 1)

julia> isequal(a + b, p / q)
true
```

While `factor(RationalPoly)` returns a `FactoredPoly` as the factorization functions, the interpretation of the results is different. First, the factors are added not multiplied. Second, the power applies only to the denominator. For example, `(-x² + 2x + 2) / (x² + x + 1) => 2` above means `(-x² + 2x + 2) / (x² + x + 1)^2`.
