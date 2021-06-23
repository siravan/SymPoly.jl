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

function test_factor(x)
    ps = [
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
        x^8 - 4x^6 + 16x^2 - 16
    ]

    k = 0
    outcome = true
    for p = ps
        try
            printstyled(p, '\n'; color=:green)
            f = factor(p)
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
        catch e
            println(e)
        end
    end
    outcome
end

#############################################################################

@testset "arith" begin
    @test test_eq(x, (p,q)->(p+q)-p-q, "add")
    @test test_eq(x, (p,q)->p-(p÷q)*q-(p%q), "mul")
    @test test_eq(x, (p,q)->p % gcd(p,q)+q % gcd(p,q), "gcd"; max_deg=5)
    @test test_deriv(x)
    @test test_factor(x)
    @test test_fraction(x)
end
