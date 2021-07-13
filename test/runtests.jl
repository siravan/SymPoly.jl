using Random

include("utils.jl")

@syms 𝑥
@polyvar x

function test_deriv(x; n=10, min_deg=1, max_deg=6, sparcity=0.5)
    k = 0
    outcome = true
    for i = 1:n
        p = generate_rand_poly(x; min_deg, max_deg, sparcity)
        printstyled("P = ", unrationalize(p), '\n'; color=:green)
        try
            q = derivative(p)
            printstyled("∂p/∂x = ", unrationalize(q), '\n'; color=:red)
            k += 1
        catch e
            println(e)
        end
    end
    outcome
end

function test_factor(x; method=:schubert_kronecker)
    ps = Any[
        x^3 + 14x^2 + 56x + 64,
        18x^3 - 57x^2 + 53x - 12,
        6x^3 - 73x^2 - 86x + 273,
        6x^4 - x^3 + 4x^2 - x - 2,
        6x^4 - 11x^3 + 8x^2 - 33x - 30,
        5x^5 - 6x^4 - 24x^3 + 20x^2 + 7x - 2,
        15x^5 - 11x^4 + 47x^3 + 27x^2 - 38x + 8,
        x^4 - 4,
        x^4 - 8x^2 - 9,
        6x^4 - 7x^3 + 5x^2 - 20x + 17,
        x^6 - 1,
        x^6 + 1,
        x^5 + x + 1,
        6x^2 + 11x - 35,
        25x^4 - 16,
        6x^4 - 19x^3 + 24x^2 - 13x + 4,
        2x^6 + 2x^5 + 2x^4 + 4x^3 + 5x^2 - 3x - 2,
        8x^5 - 48x^4 + 90x^3 - 90x^2 + 117x - 27,
        # 30x^5 + 39x^4 + 35x^3 + 25x^2 + 9x + 2,
        x^5 - x^4 - 2x^3 + 2x^2 + x - 1,
        x^10 + 2x^9 + 2x^8 + 2x^7 + x^6 + x^5 + 2x^4 + x^3 + x^2 + 2x + 1,
        x^8 - 4x^6 + 16x^2 - 16,

        𝑥^3 + 14𝑥^2 + 56𝑥 + 64,
        18𝑥^3 - 57𝑥^2 + 53𝑥 - 12,
        6𝑥^3 - 73𝑥^2 - 86𝑥 + 273,
        6𝑥^4 - 𝑥^3 + 4𝑥^2 - 𝑥 - 2,
        6𝑥^4 - 11𝑥^3 + 8𝑥^2 - 33𝑥 - 30,
        5𝑥^5 - 6𝑥^4 - 24𝑥^3 + 20𝑥^2 + 7𝑥 - 2,
        15𝑥^5 - 11𝑥^4 + 47𝑥^3 + 27𝑥^2 - 38𝑥 + 8,
        𝑥^4 - 4,
        𝑥^4 - 8𝑥^2 - 9,
        6𝑥^4 - 7𝑥^3 + 5𝑥^2 - 20𝑥 + 17,
        𝑥^6 - 1,
        𝑥^6 + 1,
        𝑥^5 + 𝑥 + 1,
        6𝑥^2 + 11𝑥 - 35,
        25𝑥^4 - 16,
        6𝑥^4 - 19𝑥^3 + 24𝑥^2 - 13𝑥 + 4,
        2𝑥^6 + 2𝑥^5 + 2𝑥^4 + 4𝑥^3 + 5𝑥^2 - 3𝑥 - 2,
        8𝑥^5 - 48𝑥^4 + 90𝑥^3 - 90𝑥^2 + 117𝑥 - 27,
        # 30𝑥^5 + 39𝑥^4 + 35𝑥^3 + 25𝑥^2 + 9𝑥 + 2,
        𝑥^5 - 𝑥^4 - 2𝑥^3 + 2𝑥^2 + 𝑥 - 1,
        𝑥^10 + 2𝑥^9 + 2𝑥^8 + 2𝑥^7 + 𝑥^6 + 𝑥^5 + 2𝑥^4 + 𝑥^3 + 𝑥^2 + 2𝑥 + 1,
        𝑥^8 - 4𝑥^6 + 16𝑥^2 - 16,
    ]

    k = 0
    outcome = true
    for p = ps
        try
            printstyled(p, '\n'; color=:green)
            f = factor(p; method=method)
            outcome = outcome && (f != nothing)
            printstyled(poly(f), '\n'; color=:red)
            println(f)
            k += 1
        catch e
            println(e)
        end
    end
    outcome
end

function test_fraction(x; n=10)
    k = 0
    outcome = true

    for i = 1:n
        p = generate_rand_poly(x; min_deg=1, max_deg=6)

        q = generate_rand_poly(x; min_deg=1, max_deg=2) *
            generate_rand_poly(x; min_deg=1, max_deg=2) *
            generate_rand_poly(x; min_deg=2, max_deg=3)

        printstyled("P/Q = ", p/q, '\n'; color=:green)

        try
            f = factor(p / q)
            print(f, '\n')
            if !iszero(poly(f) - p / q)
                outcome = false
            else
                k += 1
            end
            f = factor(sym(p, x => 𝑥), sym(q, x => 𝑥))
            printstyled(f, '\n'; color=:red)
        catch e
            println(e)
        end
    end
    outcome
end

function test_factor_extended(x, n, d=5; method=:roots_comb, seed=0)
    frags = [
        (x, i) -> (x - i),
        (x, i) -> (rand(2:11)*x - i),
        (x, i) -> (x^i + rand(1:100)),
        (x, i) -> (rand(1:20)*x^i + rand(1:100)),
        (x, i) -> (x^i - 1),
        (x, i) -> (x^(2^i) - i),
        (x, i) -> (rand(2:11)*x^(2^i) + i),
        (x, i) -> (x^2 + x + 1),
        (x, i) -> (x + 1),
        (x, i) -> (rand(1:20)*x + 1),
        (x, i) -> (x^3 + 7x^2 + 13x + 1),
        (x, i) -> (x^5 + x^4 + 11x^3 + 13x^2 + 5x + 1),
        (x, i) -> (rand(2:20)*x^3 + 7x^2 + 13x + 1),
        (x, i) -> (rand(2:20)*x^5 + x^4 + 11x^3 + 13x^2 + 5x + 1),
        (x, i) -> (rand(2:20)*x^7 + 1),
        (x, i) -> sum(rand(3:10)*x^j for j=1:rand(5:10)),
    ]

    if seed != 0
        Random.seed!(seed)
    end

    m = length(frags)

    for j = 1:n
        println("-------------------------------------")
        p = one(x)
        n₁ = 0
        for i=1:d
            q = frags[rand(1:m)](x, i)
            k = rand(1:2)
            print('(', q, ")^", k)
            if i < d
                print(" * ")
            else
                println()
            end
            p *= q^k
            n₁ += k
        end

        printstyled(p, "\n"; color=:blue)
        f = factor(p; method=method)
        printstyled(f, "\n"; color=:magenta)

        n₂ = sum(last(v) for v in f.factors)

        if p isa AbstractPolynomialLike
            println("Δ = ", poly(f) - p)
        else
            println("Δ = ", expand(simplify(f) - p))
        end

        if n₂ < n₁
            printstyled("possible incomplete factorization\n"; color=:yellow)
        else
            printstyled("OK!\n"; color=:green)
        end
    end
    true
end

function swinnerton_dyer(x, n)
    a = zeros(n)
    r = 2
    for i = 1:n
        a[i] = sqrt(r)
        r = nextprime(r+1)
    end

    p = one(x)
    for i = 0:2^n-1
        k = sum(isodd(i >> (j-1)) ? a[j] : -a[j] for j = 1:n)
        p *= x + k
    end

    return last(integer_poly(p))
end

#############################################################################

function test_all()
    @test test_eq(x, (p,q)->(p+q)-p-q, "add")
    @test test_eq(x, (p,q)->p-(p÷q)*q-(p%q), "mul")
    @test test_eq(x, (p,q)->p % gcd(p,q)+q % gcd(p,q), "gcd"; max_deg=5)
    @test test_deriv(x)
    println("********* Schubert-Kronecker ************")
    @test test_factor(x; method=:schubert_kronecker)
    println("********** Roundabout *******************")
    @test test_factor(x; method=:roundabout)
    println("********* Roots Combinations*************")
    @test test_factor(x; method=:roots_comb)
    @test test_factor_extended(x, 50; method=:roots_comb)
    @test test_factor_extended(𝑥, 50; method=:roots_comb)
    @test test_fraction(x)
end

@testset "arith" begin test_all() end
