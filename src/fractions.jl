
struct FactoredFraction
    factors::Array{Tuple{Poly,Poly}}

    FactoredFraction() = new(Array{Tuple{Poly,Poly}}[])
end

factors(f::FactoredFraction) = f.factors

function add_factor!(f::FactoredFraction, p::Poly, q::Poly)
    push!(f.factors, (p, q))
end

sym(f::FactoredFraction) = sum(map(v -> sym(first(v)) / sym(last(v)), factors(f)); init=0)
poly(f::FactoredFraction) = sum(map(v -> first(v) / last(v), factors(f)); init=0)

Base.show(io::IO, f::FactoredFraction) = print(io, sym(f))

##############################################################################

function expand_frac(p::Poly, q::Poly)
    !isequal(var(p), var(q)) && error("the numerator and denominator should have the same main variable")

    p = rationalize(p)
    q = rationalize(q)

    h = FactoredFraction()
    p₀, p = divide_poly(p, q)

    if !iszero(p₀)
        add_factor!(h, p₀, polyas(p₀, 1))
    end

    f = decompose(q)
    q₂ = copy(q)

    for w in factors(f)
        v, k = first(w), last(w)
        q₁ = v ^ k
        q₂ = q₂ / q₁
        g, s, t = gcdx(q₁, q₂)
        add_factor!(h, t*p % q₁, q₁)
        p = s * p
    end

    h
end
