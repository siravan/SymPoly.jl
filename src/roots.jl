"""
    solve_newton is a symbolic Newton-Ralphson solver
    f is a symbolic equation to be solved (f ~ 0)
    x is the variable to solve
    x₀ is the initial guess
"""
function solve_newton(f, x, x₀; abstol=1e-8, maxiter=50)
    xₙ = Complex(x₀)
    ∂f = differentiate(f, x)

    for i = 1:maxiter
        xₙ₊₁ = xₙ - f(xₙ) / ∂f(xₙ)

        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return nothing
end

"""
    find_roots returns all the real and complex roots of a polynomial (poly)
    for variable x

    the output is rs, sz, where rs is a list of real roots and zs is a list of
    complex roots

    if poly is not a polynomial, the number of desired roots (n) needs to be
    provided. For example, find_roots(sin(x), x, 10) returns 10 different, but
    not necessarily sequential, multiples of π
"""
function find_roots(p::AbstractPolynomial, x, n=-1; abstol=1e-8)
    n = (n == -1 ? deg(p, x) : n)
    rs = Float64[]
    zs = Complex[]

    while n > 0
        z = solve_newton(p, x, exp(2π*im*rand()))
        if z != nothing            
            if abs(imag(z)) < abstol
                r = real(z)
                push!(rs, r)
                p = p ÷ (x - r)
                n -= 1
            else
                if abs(real(z)) < abstol
                    z = Complex(0, imag(z))
                end
                append!(zs, [z, conj(z)])
                p = p ÷ (x^2 - 2x*real(z) + abs2(z))
                n -= 2
            end
        end
    end
    sort(rs), zs
end
