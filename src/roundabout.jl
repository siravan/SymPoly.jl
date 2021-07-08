function landau_mignotte(p::AbstractPolynomial)
    l = ceil(Int, hypot(coefficients(p)...))
    d = deg(p) ÷ 2
    maximum([binomial(d-1,j)*l + binomial(d-1,j-1)*leading(p) for j = 1:d])
end

function factor_distinct_degree(p::Polynomial{true, ℤₚ{n}}) where n
    f = FactoredPoly()
    p₀ = p
    x = modular(n, var(p))
    d = 1
    while deg(p) > 1 # && d <= deg(p) ÷ 2
        # q = modular(n, x^(n^d) - x)
        q = modpow(x, n^d, p₀) - x
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

@polyvar 𝑢

function factor_equal_degree(p::Polynomial{true, ℤₚ{n}}, d, i₀=1; maxiter=100) where n
    f = FactoredPoly()
    p₀ = p
    # x = var(p)
    x = var(p)
    e = (n^d - 1) ÷ 2

    for i = i₀:n
        a = modular(n, a+i)
        q = modpow(a, e, p₀) - one(p₀)
        g = gcd(p, q)
        # qᵤ = modular(n, 𝑢^e - 1)
        # pᵤ = p(x => 𝑢 - i)
        # gᵤ = gcd(pᵤ, qᵤ)
        # g = gᵤ(𝑢 => x + i)
        g /= leading(g)

        if !isone(g)
            if deg(g) == d
                add_factor!(f, g, 1)
            else
                f₂ = factor_equal_degree(g, d, i+1)
                combine!(f, f₂)
            end
            p = remove_factor(p, g)
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
                        p = remove_factor(p, ρ)
                        for i in a
                            push!(S, i)
                        end
                    end
                end
            end
        end
        d += 1
    end

    push!(fs, p)
    fs
end

function factor_monic(p::AbstractPolynomial, n)
    # n = nextprime(landau_mignotte(p)*2)
    P = modular(n, p)
    f = factor_distinct_degree(P)
    f = factor_equal_degree(f)
    [demodular(v) for v in first.(factors(f))]
end

# p should be integer and monic
function factor_roundabout(p::AbstractPolynomial, n)
    !isprime(n) && error("$n is not a prime!")
    λ = landau_mignotte(p)
    @info "λ = $λ"
    lc = leading(p)
    p, undo = standard_form(p)
    x = var(p)

    f = FactoredPoly()
    p₁ = one(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            f₁ = factor_monic(v, n)
            println("f₁: ", f₁)
            f₂ = lift(v, f₁, n, λ)
            # f₂ = f₁
            println("f₂: ", f₂)
            f₃ = find_integer_factorization(v, f₂)
            println("f₃: ", f₃)

            for u in f₃
                if deg(u) > 0
                    w₁ = undo(u)
                    add_factor!(f, w₁, k)
                    p₁ = p₁ * w₁^k
                end
            end
        end
    end


    ρ = lc // leading(p₁)
    if denominator(ρ) == 1
        ρ = numerator(ρ)
    end

    if !isone(ρ)
        add_factor!(f, ρ, 1)
    end

    f
end

function factor_roundabout(eq, n)
    p, v = prewrap(eq)
    unwrap(factor_roundabout(p, n), v)
end


##############################################################################

function remove_factor(p::AbstractPolynomialLike, f)
    q, r = divrem(rationalize(p), f)
    !iszero(r) && error("$f is not a proper factor of $p")
    last(integer_poly(q))
end

function lift(p::AbstractPolynomial, s₁, t₁, n, λ)
    sₙ = modular(n, s₁)
    tₙ = modular(n, t₁)
    _, σₙ, τₙ = gcdx(sₙ, tₙ)

    s, t = s₁, t₁
    d = p - s*t
    cₙ = modular(n, d ÷ n)

    i = 2
    while n^i < 2λ && !iszero(d)
        s̄ₙ = rem(τₙ*cₙ, modular(n, s))
        t̄ₙ = rem(σₙ*cₙ, modular(n, t))
        s += n^(i-1) * demodular(s̄ₙ)
        t += n^(i-1) * demodular(t̄ₙ)
        d = p - s*t
        cₘ = modular(n, d ÷ n^i)
        i += 1
    end

    f = FactoredPoly()
    add_factor!(f, s)
    add_factor!(f, t)
    f
end

function lift(p::AbstractPolynomial, f, m, λ)
    nf = length(f)
    if nf <= 1
        return f
    elseif nf == 2
        s, t = f[1], f[2]
        return lift(p, s, t, m, λ)
    else
        s, t = f[1], last(integer_poly(prod(f[i] for i=2:nf)))
        h₁ = lift(p, s, t, m, λ)
        s = h₁[1]
        h₂ = lift(remove_factor(p,s), f[2:nf], m, λ)
        add_factor!(h₂, s)
        return h₂
    end
end

function modpow(p::AbstractPolynomialLike, k::Integer, q::AbstractPolynomialLike)
    x = rem(p, q)
    r = one(p)
    while k > 0
        if isodd(k)
            r = rem(r * x, q)
        end
        x = rem(x*x, q)
        k >>= 1
    end
    r
end
