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

# function factor_equal_degree(p::Polynomial{true, ℤₚ{n}}, d, i₀=0) where n
#     f = FactoredPoly()
#     p₀ = p
#     # x = var(p)
#     x = var(p)
#     e = (n^d - 1) ÷ 2
#
#     for i = i₀:n
#         a = modular(n, x+i)
#         # a = modular(n, ((i-1)÷n)*x^2 + x + ((i-1)%n))
#         q = modpow(a, e, p₀) - one(p₀)
#         g = gcd(p, q)
#         # qᵤ = modular(n, 𝑢^e - 1)
#         # pᵤ = p(x => 𝑢 - i)
#         # gᵤ = gcd(pᵤ, qᵤ)
#         # g = gᵤ(𝑢 => x + i)
#         g ÷= leading(g)
#
#         if !isone(g)
#             if deg(g) == d
#                 add_factor!(f, g, 1)
#             else
#                 f₂ = factor_equal_degree(g, d, i+1)
#                 combine!(f, f₂)
#             end
#             p ÷= g
#         end
#
#         isone(p) && return f
#     end
#     return f
# end
#
# function factor_equal_degree(f::FactoredPoly)
#     h = FactoredPoly()
#
#     for w in factors(f)
#         p, d = first(w), last(w)
#         if deg(p) == d
#             add_factor!(h, p, 1)
#         else
#             combine!(h, factor_equal_degree(p, d))
#         end
#     end
#     h
# end

function is_square_free(p::AbstractPolynomial)
    p′ = derivative(p)
    g = gcd(p, p′)
    return isone(g)
end

# function find_integer_factorization(p::AbstractPolynomial, f)
#     m = length(f)
#     S = Set{Int}()
#     fs = []
#
#     d = 1
#     while length(S) < m && d <= deg(p)÷2
#         for a in Iterators.product([1:m for i=1:d]...)
#             if length(unique(a)) == d
#                 if !any(i ∈ S for i in a)
#                     ρ = prod(f[i] for i in a)
#                     if iszero(p % ρ)
#                         push!(fs, ρ)
#                         p = remove_factor(p, ρ)
#                         for i in a
#                             push!(S, i)
#                         end
#                     end
#                 end
#             end
#         end
#         d += 1
#     end
#
#     push!(fs, p)
#     fs
# end

function factor_modular(p::AbstractPolynomial, n)
    q = modular(n, p)
    # f = factor_distinct_degree(q)
    # f = factor_equal_degree(f)
    f = factor_modular(q)
    [demodular(v) for v in first.(factors(f))]
end

# p should be integer and monic
function factor_roundabout(p::AbstractPolynomial, n)
    !isprime(n) && error("$n is not a prime!")
    lc = leading(p)
    # p, undo = standard_form(p)
    undo = x->x
    λ = landau_mignotte(p)
    x = var(p)

    f = FactoredPoly()
    p₁ = one(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            f₁ = factor_modular(v, n)
            f₂ = lift(v, f₁, n, 2λ)
            f₃ = find_integer_factorization(v, f₂)

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

function factor_roundabout(p::AbstractPolynomial; n=3, max_prime=50)
    n >= max_prime && return nothing
    n = nextprime(n)

    f₁ = factor_roundabout(p, n)

    i = 1
    while length(f₁) == 1 && i < 4
        n = nextprime(n+1)
        f₁ = factor_roundabout(p, n)
        i += 1
    end

    length(f₁) == 1 && return f₁

    f₂ = FactoredPoly()
    k = 1

    for w₁ in factors(f₁)
        v₁, k₁ = first(w₁), last(w₁)
        if deg(v₁) == 0
            k *= v₁^k₁
        elseif deg(v₁) == 1
            k *= cont(v₁)^k₁
            add_factor!(f₂, prim(v₁), k₁)
        else
            f₃ = factor_roundabout(v₁; n=n+1, max_prime=max_prime)
            if f₃ != nothing
                for w₃ in factors(f₃)
                    v₃, k₃ = first(w₃), last(w₃)
                    if deg(v₃) == 0
                        k *= v₃^(k₁*k₃)
                    else
                        add_factor!(f₂, v₃, k₁*k₃)
                    end
                end
            end
        end
    end

    if !isone(k)
        add_factor!(f₂, k, 1)
    end

    return f₂
end

factor_roundabout(eq) = wrap(factor_roundabout, eq)

##############################################################################

function remove_factor(p::AbstractPolynomialLike, f)
    q, r = divrem(rationalize(p), f)
    !iszero(r) && error("$f is not a proper factor of $p")
    integer_poly(q)
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
        s, t = f[1], integer_poly(prod(f[i] for i=2:nf))
        h₁ = lift(p, s, t, m, λ)
        s = h₁[1]
        # h₂ = lift(remove_factor(p,s), f[2:nf], m, λ)
        h₂ = lift(t, f[2:nf], m, λ)
        add_factor!(h₂, s)
        return h₂
    end
end

function modpow(p::AbstractPolynomialLike, k::Integer, q::AbstractPolynomialLike)
    x = rem(p, q)
    r = one(p)
    while k > 0
        if isodd(k)
            r = rem(r*x, q)
        end
        x = rem(x*x, q)
        k >>= 1
    end
    r
end

##############################################################################

# p is a modular polynomial
function factor_modular(p::Polynomial{true, ℤₚ{n}}) where n
    x = modular(n, var(p))
    h = x
    v = p ÷ leading(p)
    f = FactoredPoly()
    i = 0

    while deg(v) >= 2*(i+1) && i < n
        i += 1
        h = modpow(h, n, p)
        g = gcd(h - x, v)

        if !isone(g)
            f₁ = factor_equal_degree(g, i)
            if f₁ == nothing return nothing end
            for gₛ in f₁
                if deg(gₛ) > 0
                    k = 0
                    while iszero(v % gₛ)
                        k += 1
                        v ÷= gₛ
                    end
                    if k > 0
                        add_factor!(f, gₛ, k)
                    end
                end
            end
        end
    end
    f
end

function factor_equal_degree(p::Polynomial{true, ℤₚ{n}}, d) where n
    m = deg(p)
    m == d && return [p]

    g = nothing
    i = 0
    while g == nothing || !(1 < deg(g) < m)
        g = split_equal_degree(p, d)
        i += 1
        if i == 10
            return nothing
        end
    end

    if deg(g) < m
        f₁ = factor_equal_degree(g, d)
        f₂ = factor_equal_degree(p÷g, d)
        if f₁ != nothing && f₂ != nothing
            return [f₁; f₂]
        elseif f₁ == nothing && f₂ != nothing
            return f₂
        elseif f₁ != nothing && f₂ == nothing
            return f₁
        else
            return nothing
        end
    end
    return [p]
end

function split_equal_degree(p::Polynomial{true, ℤₚ{n}}, d) where n
    x = modular(n, var(p))
    a = x + rand(0:n-1)
    g = gcd(p, a)
    !isone(g) && return g ÷ leading(g)
    b = modpow(a, (n^d-1)÷2, p)
    g = gcd(p, b - 1)
    !isone(g) && !iszero(g) && return g ÷ leading(g)
    nothing
end

function traverse_patterns(pat, mask, n, k, fun)
    bit = 1

    for i = 1:min(trailing_zeros(pat), n)
        pat₁ = pat | bit
        if (pat != pat₁) && (mask & pat₁ == 0)
            nb = count_ones(pat₁)
            if nb == k
                mask |= fun(pat₁)
            elseif nb < k
                mask = traverse_patterns(pat₁, mask, n, k, fun)
            end
        end
        bit <<= 1
    end
    return mask
end

function find_integer_factorization(p::AbstractPolynomial, f::FactoredPoly)
    p₀ = copy(p)
    x = var(p)
    Φ = first.(f.factors)
    η = last.(f.factors)
    n = length(η)
    f₁ = FactoredPoly()

    fun = function(pat)
        l = [i for i=1:n if testbit(pat, i)]
        q = demodular(prod(Φ[l]))
        k = 0
        while iszero(p % q)
            k += 1
            η[l] .-= 1
            p = integer_poly(p ÷ q)
        end
        if k > 0
            add_factor!(f₁, prim(q), k)
            return sum(1<<(i-1) for i=1:n if testbit(pat,i) && η[i]<=0; init=0)
        end
        return 0
    end

    mask = 0
    println("n = ", n)
    for k = 1:n
        # if binomial(n - count_ones(mask), k) > 100000
        #     break
        # end
        mask = traverse_patterns(0, mask, n, k, fun)
    end

    integer_poly(p₀ ÷ poly(f₁)), f₁
end

##############################################################################

function factor_combined(p::AbstractPolynomialLike; N=20, first_prime=3)
    c₀, p = integer_poly_coef(p)
    p₀ = copy(p)
    x = var(p)
    f₀ = FactoredPoly()
    fc = 0  # failure count
    a = first_prime

    while deg(p, x) > N
        println(p)
        a = nextprime(a+1)
        println(a)
        q = modular(a, p)
        # if !isone(gcd(q, derivative(q))) continue end
        println(q)
        f = factor_modular(q)
        println(f)

        if f != nothing
            p, f = find_integer_factorization(p, f)

            if length(f) == 0
                fc += 1
                if fc == 5
                    break
                end
            else
                fc = 0
            end

            combine!(f₀, f)
        end
    end

    if 1 < deg(p) <= N
        f = factor_roots_comb(p)
        combine!(f₀, f)
    end

    c, ρ = integer_poly_coef(p₀ ÷ poly(f₀))
    !isone(prim(ρ)) && add_factor!(f₀, prim(ρ), 1)
    c *= cont(ρ) * c₀

    if denominator(c) == 1
         c = numerator(c)
    end

    !isone(c) && add_factor!(f₀, c, 1)

    f₀
end
