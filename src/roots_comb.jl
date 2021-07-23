factor_roots_comb(p::AbstractPolynomialLike) = factor_main(p, factor_sqfr_roots)
factor_roots_comb(eq) = wrap(factor_roots_comb, eq)

factor_roots_SSP(p::AbstractPolynomialLike) = factor_main(p, factor_sqfr_SSP)
factor_roots_SSP(eq) = wrap(factor_roots_SSP, eq)

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

function combine_factorization_algs(p::AbstractPolynomialLike, lc, fs...; resolution=0, abstol=1e-10)
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
            if n₁*log(n₁) - k*log(k) - (n₁-k)*log(n₁-k) < 15
                mask = traverse_patterns(0, mask, n, k, cmask, test_fun)
            end
        end
    catch e
    end

    deg(p) > 0 && push!(G, p)
    F, G, lc
end

###############################################################################

near_integer(x::Float64; abstol=1e-5) = abs(x - round(x)) < abstol
near_integer(x::BigFloat; abstol=1e-10) = abs(x - round(x)) < abstol
exclude(v, l) = [v[i] for i=1:length(v) if i ∉ l]
testbit(x, n) = isodd(x >> (n-1))
expand_binary(x::Int) = [i for i=1:64-leading_zeros(x) if testbit(x,i)]
expand_binary(x::BigInt) = [i for i=1:128-leading_zeros(x) if testbit(x,i)]
expand_binary(x, n) = [i for i=1:n if testbit(x,i)]

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

frac(x) = x - floor(x)
index(x, L) = 1 + floor(Int, L * x)

struct Context{T, S}
    failures::BitMatrix
    σ::Vector{T}
    L::Int
    ls::Vector{S}
    abstol::T
    cmask::S
    cutoff::Int
end

Base.size(ctx::Context) = size(ctx.failures)
index(ctx::Context, t) = index(t, ctx.L)
add_path!(ctx::Context, path) = push!(ctx.ls, path)
@inbounds is_failure(ctx::Context, k, i) = ctx.failures[k,i]
@inbounds fail(ctx::Context, k, i) = (ctx.failures[k,i] = true)
can_grow(ctx::Context, path) = count_ones(path) + count_ones(path & ctx.cmask) < ctx.cutoff

"""
    solve_SSP is the main SSP-solving search function
    solve_SSP implements a pseudo-polynomial time dynamic programming algorithm
    the actual search phase is done by traverse_SSP
"""
function solve_SSP(σ::Array{T}, μ::Array{T}; resolution=1e-5) where T<:Real
    L = round(Int, 1 / resolution)
    σ = frac.(σ)
    m = L * round(Int, max(sum(σ),one(T)))
    n = length(σ)

    failures = falses(n, m)
    S = (n <= 64 ? Int : Int128)
    path = zero(S)
    ls = S[]

    cmask = sum(one(S)<<(i-1) for i=1:n if μ[i] != σ[i]; init=zero(S))
    cutoff = (n + count_ones(cmask)) ÷ 2

    ctx = Context{T,S}(failures, σ, L, ls, resolution, cmask, cutoff)
    traverse_SSP(ctx, 1, zero(T), path)
    return ctx.ls
end

"""
    traverse_SSP is the inner loop of the SSP-solving algorithm
    this function consumes the bulk of time of factorization and needs to be fine-tuned!
"""
function traverse_SSP(ctx::Context, k, t, path)
    n, m = size(ctx)

    if k == n + 1
        if 0 < count_ones(path) < n && near_integer(t; abstol=ctx.abstol)
            add_path!(ctx, path)
            return true
        end
        return false
    end

    i = index(ctx, t)

    if i > m return false end

    if is_failure(ctx, k, i) return false end

    if can_grow(ctx, path)
        success = traverse_SSP(ctx, k+1, t+ctx.σ[k], path | (1 << (k-1))) |
                  traverse_SSP(ctx, k+1, t, path)
    else
        success = traverse_SSP(ctx, k+1, t, path)
    end

    if !success
        fail(ctx, k, i)
    end

    success
end

##############################################################################

"""
    `factor_sqfr_SSP` factors a square-free polynomial `p` using `solve_SSP`
    `lc` is the leading coefficient
"""
function factor_sqfr_SSP(p::AbstractPolynomialLike, lc; abstol=1e-10, resolution=0)
    factor_sqfr_SSP(BigFloat, p, lc; abstol, resolution)
end

function factor_sqfr_SSP(T, p::AbstractPolynomialLike, lc; abstol=1e-10, resolution=0)
    x = var(p)
    p = rationalize(p) / leading(p)
    lc = T(1.0 / cont(p))
    p = prim(p)

    G = find_factors(T, x, p, lc; resolution)
    sort!(G; by=deg)
    G, [], lc
end

"""
    `find_factors` find candidate factors of `p`
"""
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

"""
    `find_factor` is the core of the SSP factorization algorithm
"""
function find_factor(T, x, p, lc; abstol=1e-10, resolution=0)
    if resolution == 0
        resolution = 1e-3 / sqrt(deg(p))
    end

    v = polynomial(T.(coefficients(p)), terms(p))
    r, s = find_roots(T, v, x; abstol)
    σ = T.([r; 2*real.(s[1:2:end])])
    μ = T.([r; abs2.(s[1:2:end])])
    n = length(σ)

    w = sortperm(frac.(σ); rev=true)
    σ = σ[w]
    μ = μ[w]

    ls = solve_SSP(σ*lc, μ*lc; resolution)
    ls = sort(ls; by=count_ones)

    d = zero(eltype(ls))
    F = []
    x₀ = find_test_x(p)
    α = round(BigInt, lc*p(x₀))
    # α = lc*p(im)

    for c in ls
        c = c ⊻ d
        if c > 0
            l = expand_binary(c, n)
            S = sum(σ[l]; init=0.0)

            if near_integer(S*lc)
                β = round(BigInt, lc*simple_factor(x₀, σ, μ, l))
                if !iszero(β) && !iszero(α % β) continue end
                # β = lc*simple_factor(im, σ, μ, l)
                # if !iszero(β)
                #     αβ = α * conj(β) / abs2(β)
                #     if !near_integer(real(αβ)) || !near_integer(imag(αβ))
                #         continue
                #     end
                # end

                q = simple_factor(x, σ, μ, l)
                if maximum(coefficients(q)) < typemax(Int)/10
                    q = prim(integer_poly(q*lc))
                    g = q #gcd(p, q)

                    if !isone(g) && iszero(p % g)
                        d |= c
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

"""
    `find_test_x` finds the integer with the smallest absolute value as a test number
    note that 1 is used implicitely by the core algorithm
    `find_test_x` searches integers in -1, 2, -2, 3, -3,... list and returns the
    first one that results in a nonzero `p`
"""
function find_test_x(p)
    i = 1
    while true
        x₀ = BigInt(1 + i÷2)
        x₀ = isodd(i) ? -x₀ : x₀
        if p(x₀) != 0 return x₀ end
        i += 1
    end
end
