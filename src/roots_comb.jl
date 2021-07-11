function factor_roots_comb(p::AbstractPolynomialLike)
    lc = leading(p)
    p₀ = copy(p)
    x = var(p)

    f = FactoredPoly()
    p₁ = one(x)

    for w in factors(decompose(p))
        v, k = first(w), last(w)

        if deg(v, x) > 0
            v₁ = polynomial(Float64.(coefficients(v)), terms(v))
            r, s = find_roots(v₁, x)
            s = s[1:2:end]
            f₁, d = find_order_all!(p, r, s, lc)

            length(f₁) == 0 && push!(f₁, v)

            for u in f₁
                if deg(u) > 0
                    w₁ = prim(u)
                    add_factor!(f, w₁, k)
                    p₁ = p₁ * w₁^k
                    lc ÷= d^k
                end
            end
        end
    end

    ρ = p₀ ÷ p₁

    if !isone(ρ)
        ρ = last(integer_poly(ρ))
        add_factor!(f, prim(ρ), 1)
        c = leading(p₀) // leading(poly(f))

        if denominator(c) == 1
            c = numerator(c)
        end

        if !isone(c)
            add_factor!(f, c, 1)
        end
    end

    f
end

factor_roots_comb(eq) = wrap(factor_roots_comb, eq)

near_integer(x; abstol=1e-6) = abs(x - round(x)) < abstol
exclude(v, l) = [v[i] for i=1:length(v) if i ∉ l]
testbit(x, n) = isodd(x >> (n-1))

# generates all n-bit numbers with d 1 bits
# excluding mask bits
# test each pattern by calling fun
function traverse_patterns(pat, n, d, mask, fun)
    d₁ = count_ones(pat) + 1
    bit = 1

    for i = 1:min(trailing_zeros(pat), n)
        pat₁ = pat | bit
        if (pat != pat₁) && (mask & pat₁ == 0)
            if d == d₁
                if fun(pat₁)
                    mask |= pat₁
                end
            else
                mask = traverse_patterns(pat₁, n, d, mask, fun)
            end
        end
        bit <<= 1
    end
    return mask
end

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

# function find_order_all!(p::AbstractPolynomialLike, r, s, lc)
#     x = var(p)
#     σ = [r; -2*real.(s)]
#     μ = [r; abs2.(s)]
#     n = length(μ)
#     mask = 0
#     f = AbstractPolynomialLike[]
#     d = 1
#
#     m = length(r) + 2*length(s)
#     # patterns = sort([i for i = 1:2^m-1 if count_ones(i)<=m÷2]; by=count_ones)
#     patterns = generate_test_pattern(m)
#
#     for b in patterns
#         if b & mask == 0
#             l = [i for i=1:n if testbit(b, i)]
#             S = sum(σ[l]; init=0.0)
#             M = prod(μ[l]; init=1.0)
#
#             if near_integer(S*lc) && near_integer(M*lc)
#                 q = one(x)
#                 for i in l
#                     if μ[i] == σ[i]
#                         q *= x - μ[i]
#                     else
#                         q *= x^2 + σ[i]*x + μ[i]
#                     end
#                 end
#
#                 q = last(integer_poly(q*lc))
#                 g = gcd(lc, coefficients(q)...)
#                 q = q ÷ g
#
#                 if iszero(p % q)
#                     push!(f, q)
#                     lc = g
#                     d *= g
#                     mask |= b
#                 end
#             end
#         end
#     end
#     f, d
# end

function find_order_all!(p::AbstractPolynomialLike, r, s, lc)
    x = var(p)
    σ = [r; -2*real.(s)]
    μ = [r; abs2.(s)]
    n = length(μ)
    mask = 0
    f = AbstractPolynomialLike[]
    d = 1

    m = length(r) + 2*length(s)
    # patterns = sort([i for i = 1:2^m-1 if count_ones(i)<=m÷2]; by=count_ones)
    # patterns = generate_test_pattern(m)

    test_fun = function(pat)
        l = [i for i=1:n if testbit(pat, i)]
        S = sum(σ[l]; init=0.0)
        M = prod(μ[l]; init=1.0)

        if near_integer(S*lc) && near_integer(M*lc)
            q = one(x)
            for i in l
                if μ[i] == σ[i]
                    q *= x - μ[i]
                else
                    q *= x^2 + σ[i]*x + μ[i]
                end
            end

            q = last(integer_poly(q*lc))
            g = gcd(lc, coefficients(q)...)
            q = q ÷ g

            if iszero(p % q)
                push!(f, q)
                lc = g
                d *= g
                return true
            end
        end
        false
    end

    mask = 0
    for k = 1:m÷2
        mask = traverse_patterns(0, n, k, mask, test_fun)
    end

    f, d
end
