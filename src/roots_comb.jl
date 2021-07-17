using Random

###################### Common Factoring ************************************

function factor_main(p::AbstractPolynomialLike, fs...)
    x = var(p)

    for q in factors((decompose2(p)))
        v, k = first(w), last(w)
    end
end

function combine_factorization_algs(p::AbstractPolynomialLike, f, fs...)
    F, G = f(p)

    if isempty(G) || isempty(fs)
        return F, G
    end

    G₀ = []

    for g in G
        F₁, G₁ = combine_factorization_algs(p::AbstractPolynomialLike, fs[1], fs[2:end]...)
        append!(F, F₁)
        append!(G₀, G₁)
    end

    return F, G
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

function factor_roots_comb2(p::AbstractPolynomialLike)
    c₀, p = integer_poly_coef(p)
    p₀ = copy(p)
    f = FactoredPoly()
    x = var(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            v₁ = polynomial(Float64.(coefficients(v)), terms(v))
            r, s = find_roots(v₁, var(p))
            f₁, d = find_order_all!(v, r, s[1:2:end], leading(p₀))

            for u in f₁
                if deg(u) > 0
                    add_factor!(f, prim(u), k)
                end
            end
        end
    end

    r = f[findmax(deg.(first.(factors(f))))[2]]
    if length(f) > 1 && deg(r) > 10
        f₂ = factor_roots_comb(f[length(f)])
        if length(f₂) > 1
            pop!(f.factors)
            combine!(f, f₂)
        end
    end

    c, ρ = integer_poly_coef(p₀ ÷ poly(f))
    !isone(prim(ρ)) && add_factor!(f, prim(ρ), 1)
    c *= cont(ρ) * c₀
    c = denominator(c) == 1 ? numerator(c) : c
    !isone(c) && add_factor!(f, c, 1)

    f
end

factor_roots_comb(eq) = wrap(factor_roots_comb, eq)

near_integer(x; abstol=1e-4) = abs(x - round(x)) < abstol
exclude(v, l) = [v[i] for i=1:length(v) if i ∉ l]
testbit(x, n) = isodd(x >> (n-1))

function generate_test_pattern(n)
    N = 2^(n-1) + (iseven(n) ? binomial(n,n÷2)÷2 : 0) - 1
    a = Array{Int}(undef, N)
    for i = 1:n
        a[i] = 1<<(i-1)
    end
    l = 1
    r = n
    m = r
    for k = 2:n÷2
        for i = l:r
            for j = 1:trailing_zeros(a[i])
                m += 1
                a[m] = a[i] | a[j]
            end
        end
        l = r + 1
        r = m
    end
    a
end

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

function find_order_all!(p::AbstractPolynomialLike, r, s, lc; max_deg=100000, verbose=false)
    x = var(p)
    nᵣ = length(r)
    σ = [r; 2*real.(s)]
    μ = [r; abs2.(s)]
    n = length(μ)
    f = AbstractPolynomialLike[]
    d = 1
    misses = 0

    test_fun = function(pat)
        l = [i for i=1:n if testbit(pat, i)]
        S = sum(σ[l]; init=0.0)
        M = prod(μ[l]; init=1.0)

        if near_integer(S*lc) && near_integer(M*lc)
            q = one(x)
            for i in l
                if i <= nᵣ
                    q *= x - σ[i]
                else
                    q *= x^2 - σ[i]*x + μ[i]
                end
            end

            q = integer_poly(q*lc)
            g = gcd(lc, coefficients(q)...)
            q = q ÷ g

            if iszero(p % q)
                if verbose
                    printstyled("factor = ", q, '\n'; color=:red)
                end
                push!(f, q)
                p = remove_factor(p, q)
                lc = g
                d *= g
                return true
            else
                misses += 1
                if misses == 10000
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

    if deg(p) > 0
        push!(f, p)
    end

    f, d
end


##############################################################################

function factor_lattice(p::AbstractPolynomialLike)
    c₀, p = integer_poly_coef(p)
    p₀ = copy(p)
    f = FactoredPoly()
    x = var(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            v₁ = polynomial(Float64.(coefficients(v)), terms(v))
            r, s = find_roots(v₁, var(p))
            f₁, d = find_order_all!(v, r, s[1:2:end], leading(p₀); max_deg=4)

            for u in f₁
                if deg(u) > 0
                    add_factor!(f, prim(u), k)
                end
            end
        end
    end

    u = f[findmax(deg.(first.(factors(f))))[2]] # the factor with highest degree

    h = find_candidate_factors(u, 1000)

    # if deg(u) > 4
    #     lc = leading(u)
    #     println(u)
    #     u = polynomial(Float64.(coefficients(u)), terms(u))
    #     r, s = find_roots(u, x)
    #     nᵣ = length(r)
    #     σ = [r; 2*real.(s[1:2:end])]
    #     μ = [r; abs2.(s[1:2:end])]
    #     n = length(σ)
    #     z = σ .- floor.(σ)
    #     m = sum(z)
    #
    #     k = 1.0
    #     l = collect(1:n)
    #
    #     while k < m
    #         println(m)
    #         try
    #             b, _  = subsetsum(z, k)
    #             q = one(x)
    #             for i = 1:length(b)
    #                 if b[i] > 0.5
    #                     if i <= nᵣ
    #                         q *= x - σ[l[i]]
    #                     else
    #                         q *= x^2 - σ[l[i]]*x + μ[l[i]]
    #                     end
    #                 end
    #             end
    #
    #             q = integer_poly(q*lc)
    #
    #             if iszero(u % q)
    #                 println(q)
    #                 z = z[b[1:end-1] .< 0.5]
    #                 l = l[b[1:end-1] .< 0.5]
    #             end
    #
    #         catch e
    #             println(e)
    #             k += 1.0
    #         end
    #     end
    # end

    f, h
end

function simple_factor(σ, μ, i)
    if μ[i] == σ[i]
        return x - σ[i]
    else
        return x^2 - σ[i]*x + μ[i]
    end
end

function find_candidate_factors(p::AbstractPolynomialLike, runs; max_found=1, abstol=1e-8, atol=0)
    find_candidate_factors(Float64, p, runs; max_found=max_found, abstol=abstol, atol=atol)
end

function find_candidate_factors(T, p::AbstractPolynomialLike, runs; max_found=1, abstol=1e-15, atol=0)
    p = polynomial(T.(coefficients(p)), terms(p))
    r, s = find_roots(T, p, x)
    nᵣ = length(r)
    σ = [r; 2*real.(s[1:2:end])]
    μ = [r; abs2.(s[1:2:end])]
    ζ = σ .- floor.(σ)
    n = length(σ)
    m = round(Int, sum(ζ))
    k = 1.0
    f = []
    found = 0

    for i = 1:runs
        (length(ζ) < 2 || found == max_found) && break
        θ = shuffle(1:length(ζ))
        ζ = ζ[θ]
        σ = σ[θ]
        μ = μ[θ]

        k = T(rand(1:m-1))
        b, _ = local_subsetsum(ζ, k; atol=atol)

        if !ismissing(b[1]) && b[end] != 1
            l = (b .> 0.5)
            sum_σ = sum(σ[l]; init=0.0)

            if near_integer(sum_σ) && sum(l) > 0
                q = prod(simple_factor(σ, μ, i) for i = 1:length(σ) if l[i]==1)
                q = integer_poly(q)

                println("testing ", q)

                if deg(q) > 0 && iszero(p % q)
                    if !any(isequal.(q, f))
                        push!(f, q)
                        p = remove_factor(p, q)
                        ζ = ζ[l .== 0]
                        σ = σ[l .== 0]
                        μ = μ[l .== 0]
                        m = round(Int, sum(ζ))
                        found += 1
                    end
                end
            end
        end
    end
    push!(f, integer_poly(p))
    f
end
