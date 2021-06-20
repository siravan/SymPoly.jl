# SymPoly.jl

**SymPoly.jl** is a collection of routines and algorithms to process general (arbitrary coefficients) and constant-coefficient polynomials on the top of [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl). **SymPoly** delegates most of the low-level computations to the **Symbolics** infrastructure.

## General Polynomials

### Conversion from and to Symbolics Expressions

We can convert a Symbolics.jl expression to a polynomial using `poly(eq)` function, where `eq` is a **Symbolics** expression. For examples,

```julia
  @syms x y

  eq = x^2 + 2x*y + y^2    # A Symbolics expression
  p = poly(eq, x)          # A SymPoly polynomial in variable x
```

Note that the main difference between a **Symbolics** polynomial expression and a **SymPoly** polynomial (of type `Poly` defined in `src/poly.jl`) is the presence of a main variable (`x` here). In `eq` above, `x` and `y` have the same standing; whereas `p` is a polynomial in the powers of `x` (so, the coefficient of `x^1` is `2y`).

`poly` has other forms too. If the **Symbolics** expression has only one variable, we can skip the second argument (`poly(eq)`). For example, `poly(x^3 - 3x + 7)` is valid. We can also add a coefficient type as the first argument (`poly(T, eq, x)` and `poly(T, eq)`). For example, `poly(Int, x^5 - 2x, x)` defines a polynomial with integer coefficients. An important class of polynomials are the ones with rational coefficients, which are used extensively in the factorization algorithms and have type `Rational{BigInt}`.

The conversion back from a **SymPoly** to a **Symbolics** expression is done with the help of `sym(p::Poly)` function.

```julia
  p1 = x^3 - y^2 + 3x*y + 5
  q = poly(p1, x)
  p2 = sym(q)
  println(isequal(p1, p2))  # => true
```

The following helper and accessor functions can be used to get various components of a `Poly` expression or convert it to different standard forms:

* `terms(p::Poly)`: returns a `Dict{Int, Any}` that maps power (an integer) to the corresponding coefficient (a **Symbolic** expression).
* `var(p::Poly)`: returns the main polynomial variable.
* `degree(p::Poly)`: returns the degree of a polynomial (in the main variable).
* `extract_coef(p::Poly)`: returns an array of the coefficients.
* `leading(p::Poly)`: returns the leading coefficient (the coefficient of the highest non-zero power of `x`).
* `cont(p::Poly)`: returns the *content* of a polynomial with integer or rational coefficients, i.e., the GCD of the coefficients.
* `prim(p::Poly)`: returns a *primitive* polynomial, i.e., `p / cont(p)`, so the coefficients are co-primes.

### Basic Operations

**SymPoly** supports basic operations on polynomials, including `+`, `-`, `*`, and `^`, which are equivalent to their **Symbolics** counterparts. However, division is different and performs a long division. For examples,

```julia
  @syms x

  p = x^2 - 2x + 1  # this is a Symbolics expression
  q = x - 1         # this is a Symbolics expression
  println(p/ q)       
```

returns

```julia
(1 + x^2 - (2x))*((x - 1)^-1)
```

whereas, `println(poly(p) / poly(q))` returns `x - 1.0`. Note that only one of the arguments to a binary operation needs to be of type `Poly`. For example, we can run the last expression as `poly(p) / q`.

Polynomial reminder is also defined using `%`. In addition, we can get both the quotient and reminder with one function call as `q, r = divide_poly(p, q)`.

### Greatest Common Divisor (GCD)

Finding the GCD of two polynomials is one of the basic blocks of many advanced polynomial manipulations (such as factorization). **SymPoly** provides different GCD implementations (*naive* and *monic*), which are selected based on the type of the polynomial by calling `gcd(p::Poly, q::Poly)`. For example,

```julia
  p = poly(x^2 - 2x + 1)
  q = poly(x - 1)
  g = gcd(p, q)
  println(g)
```
prints `x - 1`. We can also use `gcd` for polynomials whose coefficients are functions of other symbols:

```julia
  p = poly(x^2 - y^2, x)
  q = poly(x^2 - 2x*y + y^2, x)
  g = gcd(p, q)
  println(g)
```

The result is `2x*y - (2(y^2))`. We expect `x - y` as `p` is `(x-y)*(x+y)` and `q` is `(x-y)^2`. However, polynomial GCD has a multiplicative ambiguity. The result is equal to `2y*(x-y)`, which is a constant multiple of the expected result (note that `y` is considered constant when we define a polynomial of `x`).

The extended Euclid algorithm is available as `gcdx(p::Poly, q::Poly)` and returns three values `g, s, t`, such as `g` is the GCD and `g = s*p + t*q`.

### Polynomial Evaluation

The value of a polynomial `p` for a certain value of the main variable (`x`), say `x₀`, can be calculated as `p(x₀)`. For example,

```julia
  p = poly(3x^2 - 5x*y + 7y^2, x)
  println(p(2))
```

prints `12 + 7(y^2) - (10y)`.

### Root-Finding

We can find the real and complex roots of a polynomial with constant coefficients using `find_roots(p::Poly, x)`. It returns two arrays, the first one is a list of the real roots and the second one is a list of the complex ones. For example,

```julia
  p = (x-1)*(x-2)*(x^2+2x+5)
  r, s = find_roots(p, x)
```

returns `r = [1.0, 2.0]` and `s = Complex[-1.0 + 2.0im, -1.0 - 2.0im]`.

### Polynomial Factorization

[Polynomial factorization](https://en.wikipedia.org/wiki/Factorization_of_polynomials) is a critical operation and is required for symbolic integration, among other applications. From this point on, we are limited to polynomials with constant coefficients. Specially, we work with rational polynomials, whose coefficients are of type `Rational{BigInt}`. We can use the helper function `rationalize(p::Poly)` to convert a `Poly` expression to a rational one. Currently, three factorization algorithms are provided:

1. `factor_schubert_kronecker(p::Poly)`, which as its name implies uses the Schubert-Kronecker algorithm. Note that this algorithm is based on old ideas (Friedrich Theodor von Schubert died in 1825 and Leopold Kronecker died in 1891) and is not efficient. However, it is usable for low-degree polynomials.

2. `decompose(p::Poly)` that uses the [Yun's algorithm](https://en.wikipedia.org/wiki/Square-free_polynomial) to decompose a polynomial into a list of co-prime and square-free factors. This is a fast algorithm but the decomposition is incomplete. For some applications, such as simplifying rational expressions, this is all that is needed. In other occasions, it can serve as the first step in a completet factorization (see below). Function `recompose(f::FactoredPoly)` is the reverse function.

3. `factor(p::Poly)` is the main factorization entry point and currently uses a combination of square-free decomposition (first step) and the Schubert-Kronecker algorithm (second-step) to factor a polynomial. The plan is to switch to an efficient *Cantor–Zassenhaus algorithm* for the second step.

All three functions returns a value of type `FactoredPoly`, which is a list of `factor => power` pairs.

Let's look at an example. We start with a polynomial with known factors:

```julia
  p = expand((x-2) * (x+5)^2 * (x^2 + 2x + 5))
```

Now, `p` is `x^5 + 26(x^3) + 10(x^4) - 250 - (75x)`. We factor it as `f = factor(poly(p))` and list the factors as `factors(f)`:

```julia
3-element Vector{Pair{Poly, Int64}}:
  x - 2 => 1
  5 + x^2 + 2x => 1
  5 + x => 2  
```

We can generate a **Symbolics** expression showing the factorization as `sym(f)`, which returns `(x - 2)*(5 + x^2 + 2x)*((5 + x)^2)`. We can also convert it to a `Poly` by `poly(f)`; however, the output will be simplified to its non-factorized form.

### Partial Fraction Decomposition

[Partial Fraction Decomposition](https://en.wikipedia.org/wiki/Partial_fraction_decomposition) is the conversion of a fraction, where both the numerator and denominator are polynomials, to a sum of simpler fractions. Partial fraction decomposition has many applications, including in symbolic integration.

**SymPoly** uses the Hermite's method based on the square-free decomposition of the denominator to perform partial fraction decomposition. Let `p` be the numerator and `q` the denominator, where both are polynomials on the same variable. We can calculate the partial fraction decomposition of `p / q` by using `expand_frac(p, q)`. For example,

```julia
  p = poly(x^4 + 6x^3 + 7x^2 +6x + 4, x)
  q = poly(x^6 + 2x^5 + 4x^4 + 4x^3 + 4x^2 + 2x + 1, x)
  f = expand_frac(p, q)
```

returns `(2 + 2x - (x^2))*((1 + 2x + x^4 + 3(x^2) + 2(x^3))^-1) + 2((1 + x^2)^-1)`. Note that `f` is of type `FactoredFraction`. We can list individual factors using `factors(f)`:

```julia
2-element Vector{Tuple{Poly, Poly}}:
  (2, 1 + x^2)
  (2 + 2x - (x^2), 1 + 2x + x^4 + 3(x^2) + 2(x^3))
```

where each tuple represents a fraction, the first item is the numerator and the second one the denominator. So, `(2, 1 + x^2)` means `2 / (1 + x^2)`.
