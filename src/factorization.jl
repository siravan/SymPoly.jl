"""
    square-free decomposition of p
    Yun's algorithm

    based on Algorithm 6.2 in "Computer Algebra, Concepts and Techniques" by Edmund A. Lamangna
"""
function decompose(p::AbstractPolynomial)
    deg(p) < 1 && return FactoredPoly(false, p)
    p₀ = copy(p)
    p = rationalize(p)
    f = FactoredPoly()
    p′ = derivative(p)
    g = gcd(p, p′)
    r = p ÷ g
    i = 1
    x = var(p)

    while !isone(g) && deg(r, x) > 0
        s = gcd(g, r)
        p = r ÷ s
        if !isone(p)
            add_factor!(f, p / leading(p), i)
        end
        i += 1
        r = s
        g = g ÷ s
    end

    if !isone(r)
        add_factor!(f, r / leading(r), i)
    end

    c = p₀ ÷ poly(f)
    if !isone(c)
        add_factor!(f, c, 1)
    end

    return unrationalize(f)
end

decompose(eq) = wrap(decompose, eq)

function decompose2(p::AbstractPolynomial)
    deg(p) < 1 && return FactoredPoly(false, p)
    p₀ = copy(p)
    f = FactoredPoly()
    p′ = derivative(p)
    g = gcd(p, p′)
    r = integer_poly(p ÷ g)
    i = 1
    x = var(p)

    while !isone(g) && deg(r, x) > 0
        s = gcd(g, r)
        p = integer_poly(r ÷ s)
        if !isone(p)
            add_factor!(f, p ÷ leading(p), i)
        end
        i += 1
        r = s
        g = integer_poly(g ÷ s)
    end

    if !isone(r)
        add_factor!(f, r ÷ leading(r), i)
    end

    c = integer_poly(p₀ ÷ poly(f))
    if !isone(c)
        add_factor!(f, c, 1)
    end

    return unrationalize(f)
end

##################### Helper functions for Factorization #####################

# Lagrange interpolation algorithm for rational polynomials
function interp_poly(x, xs, ys)
    length(xs) != length(ys) && error("interp_poly needs the same number of x and y values")
    n = length(xs)
    Ps = 0
    for i = 1:n
        A = *([(x-xs[j]) for j=1:n if j!=i]...)
        B = *([(xs[i]-xs[j]) for j=1:n if j!=i]...)
        Ps += Rational{BigInt}(ys[i]) *  A / B
    end
    return Ps
end

function list_divisors(x)
    f = factor(x)
    l = [1]

    for (p, k) in f
        r = Int[]
        for i = 0:k
            append!(r, p^i .* l)
        end
        l = r
    end
    l
end

# note k is 1-based
function generate_cross_lists(D)
    n = *(length.(D)...)
    L = zeros(Int, (length(D), n))

    for k = 1:n
        j = k - 1
        for i = 1:length(D)
            n = length(D[i])
            L[i,k] = D[i][(j % n) + 1]
            j = j ÷ n
        end
    end
    return L
end

function assemble_factors(f, r, l)
    if !isone(r)
        add_factor!(f, r)
    end
    if !isone(l)
        add_factor!(f, l)
    end
    unrationalize(f)
end

########################## Main factorization algorithms #####################

"""
    Schubert-Kronecker Algorithm

    based on Algorithm 6.1 in "Computer Algebra, Concepts and Techniques" by Edmund A. Lamangna
"""
function factor_monic_sk(p::AbstractPolynomial)
    p = rationalize(p)
    r = copy(p)
    l = coef(p)
    p = integer_poly(p)
    x = var(p)
    f = FactoredPoly()
    A = Int[]
    D = Array{Int}[]

    i = 0
    while length(D) <= deg(r, x)÷2
        a = ((i+1)÷2) * (isodd(i) ? -1 : +1)
        u = convert(Int, p(a))

        if iszero(u)
            while iszero(r % (x-a))
                r = r ÷ (x-a)
                add_factor!(f, x-a)
            end
        else
            push!(A, a)
            divs = list_divisors(abs(u))
            push!(D, [divs; -divs])
        end
        i += 1
    end

    for d = 1:deg(r, x)÷2
        L = generate_cross_lists(D[1:d+1])
        xs = A[1:d+1]
        for i = 1:size(L,2)
            ys = L[:,i]
            q = interp_poly(x, xs, ys)

            while deg(q, x) > 0 && iszero(r % q)
                r = r ÷ q
                add_factor!(f, q)
                if deg(r) < 2d
                    return assemble_factors(f, r, l)
                end
            end
        end
        if deg(r, x) < 2*(d+1)
            return assemble_factors(f, r, l)
        end
    end
    assemble_factors(f, r, l)
end

function factor_schubert_kronecker(p::AbstractPolynomial)
    pᵢ = integer_poly(p)
    pₘ, c = to_monic(pᵢ)
    x = var(p)

    f = FactoredPoly()
    p₁ = one(p)

    for w in factors(decompose(pₘ))
        v, k = first(w), last(w)
        if deg(v, x) > 0
            for u in factors(factor_monic_sk(v))
                w₁ = prim(from_monic(first(u), c))
                add_factor!(f, w₁, k)
                p₁ = p₁ * w₁^k
            end
        end
    end

    q, r = divrem(p, p₁)

    !iszero(r) && @error "incorrect factorization"

    if !isone(q)
        add_factor!(f, q, 1)
    end

    return unrationalize(f)
end

factor_schubert_kronecker(eq) = wrap(factor, eq)

##############################################################################

"""
    Rational Fraction Decompositionusing the Hermite's method

    based on part of Algorithm 7.1 in "Computer Algebra, Concepts and Techniques" by Edmund A. Lamangna
"""
function factor_rational(n::AbstractPolynomialLike, d::AbstractPolynomialLike)
    !isequal(var(n), var(d)) && error("the numerator and denominator should have the same main variable")

    f = FactoredPoly(true)

    n = rationalize(n)
    d = rationalize(d)
    q, r = divrem(n, d)

    if !iszero(q)
        add_factor!(f, q)
    end

    for w in factors(decompose(d))
        v, k = rationalize(first(w)), last(w)
        d₁ = v ^ k
        d = d ÷ d₁
        g, s, t = gcdx(d₁, d)

        num = rationalize((t*r) % d₁)
        denom = v

        if !iszero(num)
            c = cont(num) // cont(denom)
            add_factor!(f, c * prim(num) / prim(denom), k)
        end
        r = s * r
    end

    f
end

factor_rational(r::RationalPoly) = factor_rational(numerator(r) / denominator(r))
factor_rational(p, q) = wrap(factor, p, q)
