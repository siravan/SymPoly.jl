mutable struct Wrapper
    p::AbstractPolynomial
    x::Union{Nothing, AbstractVariable}
end

wrap(p, x) = Wrapper(p, x)
wrap(p) = Wrapper(p, var(p))

unwrap(w::Wrapper) = w.p
unwrap(x::Number) = x

poly(w::Wrapper) = w.p
var(w::Wrapper) = w.x

Base.one(w::Wrapper) = one(w.p)
Base.zero(w::Wrapper) = zero(w.p)

leading(w::Wrapper) = leading(w.p)
cont(w::Wrapper) = cont(w.p)
prim(w::Wrapper) = prim(w.p)

(w::Wrapper)(x₀) = w.p(x₀)

Base.getindex(w::Wrapper, k::Integer) = coefficient(w.p, w.x^k)

function Base.setindex!(w::Wrapper, term, k::Number)
    t = w.x ^ k
    w.p += (term - coefficient(w.p, t)) * t
end

function Base.:+(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p + v.p, u.x)
    else
        return wrap(u.p + v.p, nothing)
    end
end

Base.:+(u::Wrapper, v) = wrap(u.p + v, u.x)
Base.:+(u, v::Wrapper) = wrap(u + v.p, v.x)

function Base.:-(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p - v.p, u.x)
    else
        return wrap(u.p - v.p, nothing)
    end
end

Base.:-(u::Wrapper, v) = wrap(u.p - v, u.x)
Base.:-(u, v::Wrapper) = wrap(u - v.p, v.x)

function Base.:*(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p * v.p, u.x)
    else
        return wrap(u.p * v.p, nothing)
    end
end

Base.:*(u::Wrapper, v) = wrap(u.p * v, u.x)
Base.:*(u, v::Wrapper) = wrap(u * v.p, v.x)

function Base.:/(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p / v.p, u.x)
    else
        return wrap(u.p / v.p, nothing)
    end
end

Base.:/(u::Wrapper, v) = wrap(u.p / v, u.x)
Base.:/(u, v::Wrapper) = wrap(u / v.p, v.x)

function Base.:÷(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p ÷ v.p, u.x)
    else
        return wrap(u.p ÷ v.p, nothing)
    end
end

Base.:÷(u::Wrapper, v) = wrap(u.p ÷ v, u.x)
Base.:÷(u, v::Wrapper) = wrap(u ÷ v.p, v.x)

function Base.:%(u::Wrapper, v::Wrapper)
    if isequal(u.x, v.x)
        return wrap(u.p % v.p, u.x)
    else
        return wrap(u.p % v.p, nothing)
    end
end

Base.:%(u::Wrapper, v) = wrap(u.p % v, u.x)
Base.:%(u, v::Wrapper) = wrap(u % v.p, v.x)

Base.:^(u::Wrapper, k) = wrap(u.p ^ k, u.x)

Base.iszero(w::Wrapper) = iszero(w.p)
deg(w::Wrapper) = maxdegree(w.p, w.x)
derivative(w::Wrapper) = derivative(w.p, w.x)

#############################################################################

function rationalize(p::AbstractPolynomial)
    polynomial(convert.(Rational{BigInt}, coefficients(p)), terms(p))
end

function unrationalize(p::AbstractPolynomial)
    t = terms(p)
    c = map(x -> x isa Rational && denominator(x)==1 ? numerator(x) : x, coefficients(p))
    polynomial(c, t)
end

function gcd_extended(u::AbstractPolynomial, v::AbstractPolynomial)
    u = wrap(rationalize(u))
    v = wrap(rationalize(v))

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
    if iszero(l)
        return unwrap(u), unwrap(sᵤ), unwrap(tᵤ)
    else
        return unwrap(u / l), unwrap(sᵤ / l), unwrap(tᵤ / l)
    end
end

Base.gcd(u::Wrapper, v::Wrapper) = wrap(gcd(u.p, v.p), u.x)
Base.gcd(u::Wrapper, v::AbstractPolynomial) = wrap(gcd(u.p, v), u.x)
Base.gcd(u::Wrapper, v) = wrap(gcd(u.p, v), u.x)
Base.gcd(u, v::Wrapper) = wrap(gcd(u, v.p), v.x)

Base.gcdx(u::AbstractPolynomial, v::AbstractPolynomial) = gcd_extended(u, v)
Base.gcdx(u::Wrapper, v::Wrapper) = wrap(gcd_extended(u.p, v.p), u.x)
Base.gcdx(u::Wrapper, v::AbstractPolynomial) = wrap(gcd_extended(u.p, v), u.x)
Base.gcdx(u::Wrapper, v) = wrap(gcd_extended(u.p, v), u.x)
Base.gcdx(u, v::Wrapper) = wrap(gcd_extended(u, v.p), v.x)

##############################################################################

leading(p::AbstractPolynomial) = leadingcoefficient(p)
cont(p::AbstractPolynomial) = gcd(coefficients(p)...) * sign(leading(p))
prim(p::AbstractPolynomial) = p / cont(p)

function var(p::AbstractPolynomial)
    vars = variables(p)
    if length(vars) == 1
        return vars[1]
    else
        error("Polynomial should have only one variable")
    end
end

derivative(p::AbstractPolynomial) = differentiate(p, var(p))
deg(p::AbstractPolynomial, x) = maxdegree(p, x)
deg(p::AbstractPolynomial) = maxdegree(p, var(p))
