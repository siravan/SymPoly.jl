include("utils.jl")

function test_deriv(x; n=10, min_deg=0, max_deg=20, sparcity=0.5)
    k = 0
    outcome = true
    for i = 1:n
        p = generate_rand_poly(x; min_deg, max_deg, sparcity)
        printstyled("P = ", p, '\n'; color=:green)
        try
            q = derivative(p)
            printstyled("∂p/∂x = ", q, '\n'; color=:red)
            Δp = p - integrate(q)
            if SymPoly.degree(update_order(Δp)) > 0
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

function test_factor(x)
    ps = [
        gen_poly(x,  64, 56, 14, 1),
        gen_poly(x,  -12, 53, -57, 18),
        gen_poly(x,  273, -86, -73, 6),
        gen_poly(x,  -2, -1, 4, -1, 6),
        gen_poly(x,  -30, -33, 8, -11, 6),
        gen_poly(x,  -2, 7, 20, -24, -6, 5),
        gen_poly(x,  8, -38, 27, 47, -11, 15),
        poly(x^4 - 4, x),
        poly(x^4 -8x^2 - 9, x),
        poly(6x^4 - 7x^3 + 5x^2 - 20x + 17, x),
        poly(x^6 - 1, x),
        poly(x^6 + 1, x),
        poly(x^5 + x + 1, x),
        gen_poly(x,  -35, 11, 6),
        gen_poly(x,  -16, 25),
        gen_poly(x,  4, -13, 24, -19, 6),
        gen_poly(x,  -2, 3, 5, 4, 2, 2, 2),
        gen_poly(x,  -27, 117, -90, 90, -48, 8),
        gen_poly(x,  2, 9, 25, 35, 39, 30),
        gen_poly(x,  -1, 1, 2, -2, -1, 1),
        gen_poly(x,  1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1),
        poly(x^8 - 4x^6 + 16x^2 - 16, x),
        # gen_poly(x,  -1, 0, 3, 2, -3, -6, 6, 3, -2, -3, 1),
    ]

    k = 0
    outcome = true
    for p = ps
        try
            printstyled(p, '\n'; color=:green)
            f = factor(p; verbose=true)
            outcome = outcome && (f != nothing)
            println('\t', f)
            k += 1
        catch e
            println(e)
        end
    end
    outcome
end

#############################################################################

@testset "arith" begin
    @test test_eq(x, (p,q)->(p+q)-p-q, "add")
    @test test_eq(x, (p,q)->p-(p/q)*q-(p%q), "mul")
    @test test_eq(x, (p,q)->p % gcd(p,q)+q % gcd(p,q), "gcd"; max_deg=5)
    @test test_deriv(x)
    @test test_factor(x)
end
