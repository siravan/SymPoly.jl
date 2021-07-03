using Primes

# based on https://discourse.julialang.org/t/arithmetic-modulo-primes/23895/3

struct ℤₚ{p} <: Number
   val::Integer

    function ℤₚ{n}(a) where {n}
        u = mod(a,n)
        new(u > n ÷ 2 ? u - n : u)
    end
end

val(x::ℤₚ{n}) where n = x.val

ℤₚ{n}(x::ℤₚ{n}) where n = x

Base.promote(x::ℤₚ{n}, y::Integer) where {n}=(x,ℤₚ{n}(y))
Base.promote(y::Integer, x::ℤₚ{n}) where {n}=(ℤₚ{n}(y),x)

Base.zero(::Type{ℤₚ{n}}) where {n} = ℤₚ{n}(0)
Base.one(::Type{ℤₚ{n}}) where {n} = ℤₚ{n}(1)

Base.:(==)(x::ℤₚ,y::ℤₚ) = (x.val==y.val)
Base.:(==)(x::ℤₚ{n}, k::Integer) where n = (val(x) == k)
Base.:(==)(k::Integer, x::ℤₚ) = (x == k)

Base.:+(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)+y.val)
Base.:*(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)*y.val)
Base.:-(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)-y.val)
Base.:-(x::ℤₚ{n}) where {n} = ℤₚ{n}(-Int(x.val))

Base.:/(x::ℤₚ{n}, y::ℤₚ{n}) where n = x * inv(y)
Base.:÷(x::ℤₚ{n}, y::ℤₚ{n}) where n = x * inv(y)
Base.:÷(x::ℤₚ{n}, y::Integer) where n = ℤₚ{n}(Int(x.val)*invmod(y,n))

Base.inv(x::ℤₚ{n}) where {n} = ℤₚ{n}(invmod(x.val,n))
Base.real(x::ℤₚ{n}) where {n} = x.val
Base.abs(x::ℤₚ{n}) where {n} = abs(x.val)
Base.gcd(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(gcd(x.val, y.val))
Base.gcd(x::ℤₚ{n}, y::Integer) where {n} = gcd(x.val, ℤₚ{n}(y))
Base.gcd(x::Integer, y::ℤₚ{n}) where {n} = gcd(ℤₚ{n}(x), y)

function Base.show(io::IO, m::ℤₚ{n}) where n
    if get(io,:limit, false)
        sub = Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
        print(io, m.val, map(x->sub[x],repr(n)))
   else
        print(io,"ℤₚ{$n}($(m.val))")
   end
end

###############################################################################

modular(n, v::AbstractArray) = [ℤₚ{n}(x) for x in v]

function modular(n::Integer, p::AbstractPolynomial)
    polynomial(modular(n, coefficients(p)), terms(p))
end

function demodular(p::Polynomial{true, ℤₚ{n}}) where n
    polynomial(val.(coefficients(p)), terms(p))
end

function info(p::Polynomial{true, ℤₚ{n}}) where n
    println(n)
end

function landau_mignotte(p::AbstractPolynomial)
    l = ceil(Int, hypot(coefficients(p)...))
    d = deg(p) ÷ 2
    maximum([binomial(d-1,j)*l + binomial(d-1,j-1)*leading(p) for j = 1:d-1])
end

function factor_distinct_degree(p::Polynomial{true, ℤₚ{n}}) where n
    f = FactoredPoly()
    x = var(p)
    d = 1
    while deg(p) > 1 # && d <= deg(p) ÷ 2
        println(p)
        q = modular(n, x^(n^d) - x)
        println(q)
        g = gcd(p, q)

        if !isone(g)
            add_factor!(f, g, d)
            p ÷= g
        end
        d += 1
    end
    # add_factor!(f, p, deg(p))
    f
end

function factor_equal_degree(p::Polynomial{true, ℤₚ{n}}, d, i₀=1; maxiter=100) where n
    f = FactoredPoly()
    x = var(p)
    e = (n^d - 1) ÷ 2

    for i = i₀:maxiter
        q = modular(n, (i ÷ n)*x^2 + x + (i % n))
        q = q^e - 1
        g = gcd(p, q)
        g /= leading(g)

        if !isone(g)
            if deg(g) == d
                add_factor!(f, g, 1)
            else
                f₂ = factor_equal_degree(g, d, i+1)
                combine!(f, f₂)
            end
            p ÷= g
        end

        isone(p) && return f
    end
    add_factor!(f, p, 1)
    return f
end

@polyvar 𝑢

function factor_equal_degree2(p::Polynomial{true, ℤₚ{n}}, d, i₀=1; maxiter=100) where n
    f = FactoredPoly()
    x = var(p)
    e = (n^d - 1) ÷ 2

    for i = i₀:n
        qᵤ = modular(n, 𝑢^e - 1)
        pᵤ = p(x => 𝑢 - i)
        gᵤ = gcd(pᵤ, qᵤ)
        g = gᵤ(𝑢 => x + i)
        g /= leading(g)

        # println(p, "\t", q, "\t", g)

        if !isone(g)
            if deg(g) == d
                add_factor!(f, g, 1)
            else
                f₂ = factor_equal_degree2(g, d, i+1)
                combine!(f, f₂)
            end
            p ÷= g
        end

        isone(p) && return f
    end
    add_factor!(f, p, 1)
    return f
end

function factor_equal_degree(f::FactoredPoly; maxiter=100)
    h = FactoredPoly()

    for w in factors(f)
        p, d = first(w), last(w)
        if deg(p) == d
            add_factor!(h, p, 1)
        else
            # combine!(h, factor_equal_degree(p, d; maxiter=maxiter))
            combine!(h, factor_equal_degree2(p, d; maxiter=maxiter))
        end
    end
    h
end

function is_square_free(p::AbstractPolynomial)
    p′ = derivative(p)
    g = gcd(p, p′)
    return isone(g)
end

function find_integer_factorization(p::AbstractPolynomial, f)
    m = length(f)
    S = Set{Int}()
    fs = []

    d = 1
    while length(S) < m && d < 3
        for a in Iterators.product([1:m for i=1:d]...)
            if length(unique(a)) == d
                if !any(i ∈ S for i in a)
                    ρ = prod(f[i] for i in a)
                    if iszero(p % ρ)
                        push!(fs, ρ)
                        p ÷= ρ
                        for i in a
                            push!(S, i)
                        end
                    end
                end
            end
        end
        d += 1
    end
    _, p = integer_poly(p)
    push!(fs, p)
    fs
end

function factor_decomposed(p::AbstractPolynomial, n)
    # n = nextprime(landau_mignotte(p)*2)
    P = modular(n, p)
    f = factor_distinct_degree(P)
    f = factor_equal_degree(f)
    f = [demodular(v) for v in first.(factors(f))]

    find_integer_factorization(p, f)
end

# p should be integer and monic
function factor_roundabout(p::AbstractPolynomial, n)
    S = standardize(p)
    p = poly(S)
    x = var(p)

    f = FactoredPoly()
    p₁ = one(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)
        if deg(v, x) > 0
            printstyled(v, '\n'; color=:red)
            D = factor_decomposed(v, n)
            for u in D
                if deg(u) > 0
                    # _, w₁ = integer_poly(prim(u))
                    w₁ = from_monic(S, u)
                    add_factor!(f, w₁, k)
                    p₁ = p₁ * w₁^k
                end
            end
        end
    end


    ρ = lc(S) // leading(p₁)

    if !isone(ρ)
        add_factor!(f, ρ, 1)
    end

    # return unrationalize(f)
    f
end

##############################################################################

mutable struct StandardForm
    poly::AbstractPolynomial
    sym
    xp
    xs
    coef
    lc
end

poly(S::StandardForm) = S.poly
var(S::StandardForm) = S.xp
coef(S::StandardForm) = S.coef
lc((S::StandardForm)) = S.lc

function standardize(p::Polynomial{true,T}) where T<:Integer
    μ, lc = to_monic(p)
    q = polynomial(numerator.(coefficients(μ)), terms(μ))
    return StandardForm(q, nothing, var(q), nothing, one(T), numerator(lc))
end

function standardize(p::Polynomial{true,T}) where T<:Rational
    coef, q = integer_poly(p)
    S = standardize(q)
    S.coef = coef
    return S
end

function standardize(p::Polynomial{true,T}) where T
    return standardize(rationalize(p))
end

function from_monic(S::StandardForm, p)
    n = deg(p)
    c = convert(Int, lc(S))
    x = var(S)
    q = p(x => c*x) * (1 // c^(n-1))
    polynomial(numerator.(coefficients(q)) .÷ cont(q), terms(q))
end

function original(S::StandardForm, p)
    return coef(S) * p
end
