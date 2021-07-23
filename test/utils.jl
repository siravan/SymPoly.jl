using SymPoly
using SymbolicUtils, Test
using DynamicPolynomials

function generate_rand_poly(x; min_deg=0, max_deg=5, sparcity=0.5)
    n = rand(min_deg:max_deg)
    p = one(x)
    r = rand(1:10) // rand(1:10)
    while true
        for i = 0:n
            if rand() > sparcity || i == n
                s = rand() < 0.5 ? +1 : -1
                p += x^i * s * r * rand(0:10)//rand(1:5)
            end
        end
        !isequal(p, 0) && return p
        return p
    end
end

function test_eq(x, f, name; n=10, min_deg=0, max_deg=10, sparcity=0.5)
    k = 0
    outcome = true
    for i = 1:n
        p = generate_rand_poly(x; max_deg, sparcity)
        q = generate_rand_poly(x; min_deg, max_deg, sparcity)

        if iszero(p) || iszero(q) continue end

        printstyled("P1 = ", p, '\n'; color=:green)
        printstyled("P2 = ", q, '\n'; color=:blue)

        try
            r = f(p, q)
            printstyled("Δ($name) = ", r, '\n'; color=:red)
            if !iszero(r)
                @warn "+/- mismatch: $r"
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
