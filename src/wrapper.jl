
function var(p::AbstractPolynomial)
    vars = variables(p)
    if length(vars) == 1
        return vars[1]
    elseif length(vars) == 0
        return nothing
    else
        error("Polynomial should have only one variable")
    end
end

function var(p::AbstractMonomial)
    vars = variables(p)
    if length(vars) == 1
        return vars[1]
    elseif length(vars) == 0
        return nothing
    else
        error("Polynomial should have only one variable")
    end
end

var(p::AbstractTerm) = variable(p)

function var(p)
    vars = get_variables(p)
    if length(vars) == 1
        return vars[1]
    elseif length(vars) == 0
        return nothing
    else
        error("Polynomial should have only one variable")
    end
end

wrap(fun, p::AbstractPolynomialLike) = fun(p)
wrap(fun, p::AbstractPolynomialLike, q::AbstractPolynomialLike) = fun(p, q)

function prewrap(eq)
    x = var(eq)
    x == nothing && return fun(eq * one(𝑦))
    p = poly(eq, x => 𝑦)
    p, 𝑦 => x
end

function wrap(fun, eq)
    p, v = prewrap(eq)
    q = fun(p)
    unwrap(q, v)
end

function wrap(fun, eq₁, eq₂)
    x₁ = var(eq₁)
    x₂ = var(eq₂)

    if x₁ == nothing && x₂ == nothing
        return fun(eq₁ * one(𝑦), eq₂ * one(𝑦))
    elseif x₁ == nothing
        p = eq₁ * one(𝑦)
        q = poly(eq₂, x₂ => 𝑦)
        x = x₂
    elseif x₂ == nothing
        p = poly(eq₁, x₁ => 𝑦)
        q = eq₂ * one(𝑦)
        x = x₁
    else
        !isequal(x₁, x₂) && error("incompatible main variables")
        p = poly(eq₁, x₁ => 𝑦)
        q = poly(eq₂, x₂ => 𝑦)
        x = x₁
    end

    unwrap(fun(p, q), 𝑦 => x)
end

unwrap(q::AbstractPolynomial, v) = sym(q, v)
unwrap(q::RationalPoly, v) = sym(numerator(q), 𝑦 => x) / sym(denominator(q), 𝑦 => x)
unwrap(q::FactoredPoly, v) = sym(q, v)
unwrap(q::Tuple, v) = map(x->unwrap(x,v), q)
unwrap(q, v) = q

###############################################################################

rationalize(p) = wrap(x->polynomial(convert.(Rational{BigInt}, coefficients(x)), terms(x)), p)

function unrationalize(p::AbstractPolynomial)
    t = terms(p)
    c = map(x -> x isa Rational && denominator(x)==1 ? numerator(x) : x, coefficients(p))
    polynomial(c, t)
end

unrationalize(p) = wrap(unrationalize, p)

function gcd_extended(u::AbstractPolynomialLike, v::AbstractPolynomialLike)
    # u = rationalize(u)
    # v = rationalize(v)

    !isequal(var(u), var(v)) && error("incompatible main variable")

    if deg(u) == 0 && deg(v) == 0
        return gcdx(cont(u), cont(v))
    end

    sᵤ = 1
    tᵤ = 0
    sᵥ = 0
    tᵥ = 1

    while !iszero(v)
        q, r = divrem(u, v)
        s, t = sᵤ - q*sᵥ, tᵤ - q*tᵥ
        u, sᵤ, tᵤ = v, sᵥ, tᵥ
        v, sᵥ, tᵥ = r, s, t
    end

    l = leading(u)
    return u / l, sᵤ / l, tᵤ / l
end

Base.gcdx(u, v) = wrap(gcd_extended, u, v)

##############################################################################

leading(p::AbstractPolynomialLike) = leadingcoefficient(p)
cont(p::AbstractPolynomialLike) = gcd(coefficients(p)...) * sign(leading(p))
prim(p::AbstractPolynomialLike) = polynomial(coefficients(p) .÷ cont(p), terms(p))
derivative(p::AbstractPolynomialLike) = differentiate(p, var(p))
deg(p::AbstractPolynomialLike, x) = maxdegree(p, x)
deg(p::AbstractPolynomialLike) = maxdegree(p, var(p))
coef(p::AbstractPolynomialLike) = lcm(denominator.(coefficients(p))...)

# leading(p::AbstractTerm) = leadingcoefficient(p)
# cont(p::AbstractTerm) = gcd(coefficients(p)...) * sign(leading(p))
# prim(p::AbstractTerm) = p ÷ cont(p)

leading(p) = wrap(leading, p)
cont(p) = wrap(cont, p)
prim(p) = wrap(prim, p)
derivative(p) = wrap(derivative, p)
deg(p) = wrap(deg, p)
coef(p) = wrap(coef, p)

function Base.:÷(p::AbstractPolynomialLike, k)
    polynomial(coefficients(p) .÷ k, terms(p))
end

###############################################################################

function integer_poly(p::Polynomial{true,T}) where T<:Rational
    l = coef(p)
    1//l*one(p), polynomial(Int.(numerator.(coefficients(l*p))), terms(p))
end

integer_poly(p::Polynomial{true,T}) where T<:Integer = (1, p)

function to_monic(p::Polynomial{true,T}) where T<:Rational
    c = leading(p)
    n = deg(p)
    d = degree.(terms(p))
    polynomial(coefficients(p) .* c.^((n-1).-d), terms(p)), c
end

function to_monic(p::Polynomial{true,T}) where T<:Integer
    q, c = to_monic(rationalize(p))
    last(integer_poly(q)), Int(numerator(c))
end

function from_monic(p::Polynomial{true,T}, c) where T<:Rational
    n = deg(p)
    d = degree.(terms(p))
    polynomial(coefficients(p) .* (1 .// c.^((n-1).-d)), terms(p))
end

function from_monic(p::Polynomial{true,T}, c) where T<:Integer
    q = from_monic(rationalize(p), Rational(c))
    last(integer_poly(q))
end

from_monic(p, c) = p

###############################################################################

function standard_form(p::Polynomial{true,T}) where T<:Integer
    q, lc = to_monic(prim(p))
    undo = s -> from_monic(s, lc)
    q, undo
end

function standard_form(p::Polynomial{true,T}) where T<:Rational
    _, q = integer_poly(p)
    standard_form(q)
end

function standard_form(eq)
    p, v = prewrap(eq)
    q, undo = standard_form(p)
    q, x -> unwrap(undo(x), v)
end
