using Random

###################### Common Factoring ************************************

function factor_main(p::AbstractPolynomialLike, fs...; resolution=1e-5, abstol=1e-10)
    c₀, p = integer_poly_coef(p)
    x = var(p)

    P = FactoredPoly()
    lc = leading(p)

    for w in factors((decompose(p)))
        v, k = first(w), last(w)
        F, G, lc = combine_factorization_algs(v, lc, fs...; abstol, resolution)

        for q in [F; G]
            add_factor!(P, prim(q), k)
        end
    end

    add_coefficient!(P, p, c₀)

    P
end

function add_coefficient!(P, p::AbstractPolynomialLike, c₀)
    c, ρ = integer_poly_coef(p ÷ poly(P))
    !isone(prim(ρ)) && add_factor!(P, prim(ρ), 1)
    c *= cont(ρ) * c₀
    c = denominator(c) == 1 ? numerator(c) : c
    !isone(c) && add_factor!(P, c, 1)
end

function combine_factorization_algs(p::AbstractPolynomialLike, lc, fs...; resolution=1e-5, abstol=1e-10)
    F₁, G₁, lc = fs[1](p, lc; resolution, abstol)

    if isempty(G₁) || length(fs) == 1
        return F₁, G₁, lc
    end

    # fs = (length(G₁) == 1 ? fs[2:end] : fs)
    fs = fs[2:end]
    G = []
    F = F₁

    for g in G₁
        F₂, G₂, lc = combine_factorization_algs(g, lc, fs...; resolution, abstol)
        append!(F, F₂)
        append!(G, G₂)
    end

    return F, G, lc
end

###############################################################################

function factor_sqfr_roots(p::AbstractPolynomialLike, lc; max_deg=100000, max_misses=10000, resolution=1e-5, abstol=1e-10)
    x = var(p)
    v = polynomial(Float64.(coefficients(p)), terms(p))
    r, s = find_roots(v, x; abstol)
    s = s[1:2:end]
    σ = [r; 2*real.(s)]
    μ = [r; abs2.(s)]
    n = length(μ)

    d = 1
    misses = 0
    F = []
    G = []

    test_fun = function(pat)
        l = [i for i=1:n if testbit(pat, i)]
        S = sum(σ[l]; init=0.0)
        M = prod(μ[l]; init=1.0)

        if near_integer(S*lc) && near_integer(M*lc)
            q = simple_factor(x, σ, μ, l)
            q = integer_poly(q*lc)
            g = gcd(lc, coefficients(q)...)
            q = q ÷ g

            if iszero(p % q)
                push!(F, q)
                p = remove_factor(p, q)
                lc = g
                d *= g
                return true
            else
                misses += 1
                if misses == max_misses
                    error("too many misses!")
                end
            end
        end
        false
    end

    m = length(r) + 2*length(s)
    mask = 0
    cmask = 2^length(μ) - 2^length(r)

    try
        for k = 1:min(m÷2, max_deg)
            n₁ = n - count_ones(mask)
            if n₁*log(n₁) - k*log(k) - (n₁-k)*log(n₁-k) < 18
                mask = traverse_patterns(0, mask, n, k, cmask, test_fun)
            end
        end
    catch e
    end

    deg(p) > 0 && push!(G, p)
    F, G, lc
end

###############################################################################

function factor_roots_comb(p::AbstractPolynomialLike)
    x = var(p)
    t = terms(p)
    c = coefficients(p)
    d = degree.(t)
    nz = sum(d .> 0)
    g = gcd(d...)

    f = FactoredPoly()

    if g > 1 && nz > 1
        q = sum(c[i]*x^(d[i] ÷ g) for i=1:length(d))
        f₁ = factor_roots_comb2(q)

        for w₁ in factors(f₁)
            v₁, k₁ = first(w₁), last(w₁)
            q = v₁(x => x^g)
            f₂ = factor_roots_comb2(q)

            for w₂ in factors(f₂)
                v₂, k₂ = first(w₂), last(w₂)
                add_factor!(f, v₂, k₁*k₂)
            end
        end
        return f
    else
        return factor_roots_comb2(p)
    end
end

factor_roots_comb(eq) = wrap(factor_roots_comb, eq)

near_integer(x::Float64; abstol=1e-5) = abs(x - round(x)) < abstol
near_integer(x::BigFloat; abstol=1e-5) = abs(x - round(x)) < abstol
exclude(v, l) = [v[i] for i=1:length(v) if i ∉ l]
testbit(x, n) = isodd(x >> (n-1))


# generates all n-bit numbers with d 1 bits
# excluding mask bits
# test each pattern by calling fun
function traverse_patterns(pat, mask, n, k, cmask, fun)
    bit = 1

    for i = 1:min(trailing_zeros(pat), n)
        pat₁ = pat | bit
        if (pat != pat₁) && (mask & pat₁ == 0)
            nb = count_ones(pat₁) + count_ones(pat₁ & cmask)
            if nb == k
                if fun(pat₁)
                    mask |= pat₁
                end
            elseif nb < k
                mask = traverse_patterns(pat₁, mask, n, k, cmask, fun)
            end
        end
        bit <<= 1
    end
    return mask
end

function simple_factor(x, σ, μ, l)
    prod(μ[i]==σ[i] ? (x - σ[i]) : x^2 - σ[i]*x + μ[i] for i in l; init=one(x))
end

simple_factor(x, σ, μ, i::Integer) = simple_factor(x, σ, μ, [i])

##############################################################################

function factor_sqfr_lattice(p::AbstractPolynomialLike, lc; max_found=1, abstol=1e-10, atol=1e-6, runs=500)
    factor_sqfr_lattice(Float64, p, lc; runs=runs, max_found=max_found, abstol=abstol, atol=atol)
end

function factor_sqfr_lattice(T, p::AbstractPolynomialLike, lc; max_found=1, abstol=1e-10, atol=0, runs=5000)
    p = rationalize(p) / leading(p)
    c = T(1.0 / cont(p))
    p = prim(p)
    v = polynomial(T.(coefficients(p)), terms(p))
    r, s = find_roots(T, v, x)
    nᵣ = length(r)
    σ = T.([r; 2*real.(s[1:2:end])])
    μ = T.([r; abs2.(s[1:2:end])])
    ζ = σ*c .- floor.(σ*c)
    n = length(σ)
    m = round(Int, sum(ζ))

    k = 1.0
    G = []
    found = 0

    for i = 1:runsp
        for k = 1:m-1
            (length(ζ) < 2 || found == max_found) && break

            b, _ = local_subsetsum(T.(ζ), T(k); atol=atol)

            if !ismissing(b[1])
                l = [i for i=1:length(σ) if b[i] .> 0.5]
                sum_σ = sum(σ[l]; init=0.0)

                if near_integer(sum_σ*c) && sum(l) > 0
                    q = simple_factor(x, σ, μ, l)
                    maximum(coefficients(q)) > typemax(Int)/10 && continue
                    q = prim(integer_poly(q*c))

                    if deg(q) > 0 && iszero(p % q)
                        if !any(isequal.(q, G))
                            push!(G, q)
                            p = remove_factor(p, q)
                            l₀ = [i for i=1:length(σ) if i ∉ l]
                            ζ = ζ[l₀]
                            σ = σ[l₀]
                            μ = μ[l₀]
                            m = round(Int, sum(ζ))
                            found += 1
                        end
                    end
                end
            end
        end

        θ = shuffle(1:length(ζ))
        ζ = ζ[θ]
        σ = σ[θ]
        μ = μ[θ]
    end

    deg(p) > 0 && push!(G, p)
    [], G, lc
end

###############################################################################

frac(x) = x - floor(x)
index(x, L) = 1 + floor(Int, L * x)

function solve_SSP(σ::Array{T}; resolution=1e-5) where T<:Real
    L = round(Int, 1 / resolution)
    σ = frac.(σ)
    m = L * round(Int, max(sum(σ),one(T)))
    n = length(σ)
    failures = falses(n, m)

    ls = Int[]
    path = (n <= 64 ? zero(Int) : zero(BigInt))

    traverse_SSP(1, zero(T), failures, σ, L, path, ls)
    return ls
end

function traverse_SSP(k, t, failures, σ, L, path, ls)
    n, m = size(failures)
    i = index(t, L)

    if i > m
        return false
    end

    if k == n + 1
        if 0 < count_ones(path) < n && near_integer(t; abstol=1.0/L)
            push!(ls, path)
            return true
        end
        return false
    end

    if failures[k,i]
        return false
    end

    success = traverse_SSP(k+1, t+σ[k], failures, σ, L, path | (1 << (k-1)), ls) |
              traverse_SSP(k+1, t, failures, σ, L, path, ls)

    if !success
        failures[k,i] = true
    end

    success
end

expand_binary(x, n) = [i for i=1:n if testbit(x,i)]

function factor_sqfr_SSP(p::AbstractPolynomialLike, lc; abstol=1e-10, resolution=1e-5)
    factor_sqfr_SSP(BigFloat, p, lc; abstol, resolution)
end

function factor_sqfr_SSP(T, p::AbstractPolynomialLike, lc; abstol=1e-10, resolution=1e-5)
    x = var(p)
    p = rationalize(p) / leading(p)
    lc = T(1.0 / cont(p))
    p = prim(p)

    G = find_factors(T, x, p, lc; resolution)
    sort!(G; by=deg)
    [], G, lc
end

function find_factors(T, x, p, lc; resolution=1e-5)
    F = find_factor(T, x, p, lc; resolution)
    length(F) == 1 && return F

    G = []
    for q in F
        Q = find_factors(T, x, q, lc; resolution)
        append!(G, Q)
    end
    return G
end

function find_factor(T, x, p, lc; abstol=1e-10, resolution=1e-5)
    v = polynomial(T.(coefficients(p)), terms(p))
    r, s = find_roots(T, v, x; abstol)
    σ = T.([r; 2*real.(s[1:2:end])])
    μ = T.([r; abs2.(s[1:2:end])])
    n = length(σ)

    w = sortperm(frac.(σ); rev=true)
    σ = σ[w]
    μ = μ[w]

    ls = solve_SSP(σ*lc; resolution)
    F = []

    for c in ls
        if c > 0
            l = expand_binary(c, n)
            S = sum(σ[l]; init=0.0)

            if near_integer(S*lc)
                q = simple_factor(x, σ, μ, l)
                if maximum(coefficients(q)) < typemax(Int)/10
                    q = prim(integer_poly(q*lc))
                    g = gcd(p, q)

                    if !isone(g) && iszero(p % g)
                        p = prim(remove_factor(p, g))
                        push!(F, g)
                    end
                end
            end
        end
    end

    !isone(p) && push!(F, p)
    return F
end
