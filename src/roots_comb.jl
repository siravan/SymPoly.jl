
function factor_roots_comb(p::AbstractPolynomialLike; run=1)    
    c₀, p = integer_poly(p)
    p₀ = copy(p)
    f = FactoredPoly()
    x = var(p)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            v₁ = polynomial(Float64.(coefficients(v)), terms(v))
            r, s = find_roots(v₁, var(p))
            f₁, d = find_order_all!(v, r, s[1:2:end], leading(p₀))

            # length(f₁) == 0 && push!(f₁, v)

            for u in f₁
                if deg(u) > 0
                    add_factor!(f, prim(u), k)
                end
            end
        end
    end

    c, ρ = integer_poly(p₀ ÷ poly(f))

    if deg(ρ) > 0 && run < 3
        f₂ = factor_roots_comb(ρ; run=run+1)
        combine!(f, f₂)
        c, ρ = integer_poly(p₀ ÷ poly(f))
    end

    !isone(prim(ρ)) && add_factor!(f, prim(ρ), 1)
    c *= cont(ρ) * c₀

    if denominator(c) == 1
         c = numerator(c)
    end

    !isone(c) && add_factor!(f, c, 1)

    f
end

factor_roots_comb(eq) = wrap(factor_roots_comb, eq)

near_integer(x; abstol=1e-3) = abs(x - round(x)) < abstol
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

function find_order_all!(p::AbstractPolynomialLike, r, s, lc; verbose=false)
    x = var(p)
    nᵣ = length(r)
    σ = [r; 2*real.(s)]
    μ = [r; abs2.(s)]
    n = length(μ)
    f = AbstractPolynomialLike[]
    d = 1

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

            q = last(integer_poly(q*lc))
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
            end
        end
        false
    end

    m = length(r) + 2*length(s)
    mask = 0
    cmask = 2^length(μ) - 2^length(r)

    for k = 1:m÷2
         mask = traverse_patterns(0, mask, n, k, cmask, test_fun)
    end

    if deg(p) > 0
        push!(f, p)
    end

    f, d
end
