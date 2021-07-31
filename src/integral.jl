using LinearAlgebra
using SpecialFunctions

isdependent(eq, x) = !isequal(expand_derivatives(Differential(x)(eq)), 0)

kernelize(eq, x) = x
kernelize(eq::Num, x) = kernelize(value(eq), x)
kernelize(eq::T, x) where T<:Number = one(x)
kernelize(eq::SymbolicUtils.Term, x) = (isdependent(eq,x) ? eq : one(x))
kernelize(eq::SymbolicUtils.Pow, x) = (isdependent(eq,x) && arguments(eq)[2] < 0 ? eq*log(arguments(eq)[1]) : eq)
kernelize(eq::SymbolicUtils.Mul, x) = prod(kernelize(t,x) for t in arguments(eq); init=one(x))
kernelize(eq::SymbolicUtils.Add, x) = sum(kernelize(t,x) for t in arguments(eq) if isdependent(t,x); init=zero(x))

function kernelize(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) return eq end

    p = arguments(eq)[1]
    k = arguments(eq)[2]

    if k >= 0
        return eq * p
    else
        if !is_poly(p, x) return eq*log(p) end
        q = poly(p)
        r, s = find_roots(q, var(q))
        r = nice_parameters(r)
        s = Complex.(nice_parameters(real.(s[1:2:end])), nice_parameters(imag(s[1:2:end])))

        q = sum(log(x - u) for u in r; init=zero(x)) +
            sum(atan((x - real(u))/imag(u)) for u in s; init=zero(x)) +
            sum(log(x^2 - 2*real(u)*x + imag(u)) for u in s; init=zero(x))

        return eq * p * q
    end
end

candidates(eq, x) = isdependent(eq,x) ? [eq] : []
candidates(eq::Num, x) = candidates(value(eq), x)
candidates(eq::SymbolicUtils.Add, x) = unique(∪([candidates(t,x) for t in arguments(eq)]...))

function candidates(eq::SymbolicUtils.Mul, x)
    l = Any[one(x)]
    terms = [candidates(q,x) for q in arguments(eq)]
    n = length(terms)

    for j = 1:n
        m = length(l)
        for t in terms[j]
            for k = 1:m
                push!(l, l[k]*t)
            end
        end
    end

    unique(l[2:end])
end

function candidates(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) return [one(x)] end
    p = arguments(eq)[1]
    k = arguments(eq)[2]

    # if k > 0
    #     return [p^j for j=k+1:-1:0]
    # else
    #     return [j≈-1 ? log(p) : p^j for j=k-1:0]
    # end

    if k == -1
        return [p^-1, p, log(p)]
    elseif !(k ≈ round(k)) && (2k ≈ round(2k))
        return [p^j for j=k-1:k+2]
    elseif k > 1
        # return [p^j for j=1:k+1]
        return [p^j for j=k+1:-1:0]
    else
        return [p^k, p^(k+1)]
    end
end

function integrate(eq; abstol=1e-3, num_trials=5, num_zeros=5, lo=-5.0, hi=5.0, show_basis=false)
    x = var(eq)
    if x == nothing
        @syms 𝑥
        return 𝑥 * eq
    end
    integrate(eq, x; abstol, num_trials, num_zeros, lo, hi, show_basis)
end

function integrate(eq, x; abstol=1e-5, num_trials=5, num_zeros=5, lo=-5.0, hi=5.0, show_basis=false)
    D = Differential(x)

    for i = 1:num_trials
        y = try_integrate(eq, x, i; abstol, lo, hi, show_basis)
        h = expand_derivatives(D(y)) - eq
        j = 1
        k = 0
        while j <= num_zeros
            x₀ = rand()*(hi - lo) + lo
            try
                ϵ = substitute(h, Dict(x => x₀))
                if abs(ϵ) < abstol
                    k +=1
                    if k == num_zeros
                        return y
                    end
                end
                j += 1
            catch e
                # println(e)
            end
        end
        abstol *= 0.1
    end
    @warn "no solution is found"
    nothing
end

function try_integrate(eq, x, k; abstol=1e-5, lo=-5.0, hi=5.0, show_basis=true)
    eq₁ = apply_integration_rules(eq)
    basis = generate_basis(eq₁, x, k)
    if show_basis println(basis) end
    D = Differential(x)
    Δbasis = [expand_derivatives(D(f)) for f in basis]

    n = length(basis)
    # A is an nxn matrix holding the values of the fragments at n random points
    # b hold the value of the input function at those points
    A = zeros(n, n)
    b = zeros(n)

    i = 1
    while i <= n
        x₀ = rand()*(hi - lo) + lo
        d = Dict(x => x₀)
        try
            for j = 1:n
                A[i, j] = Float64(substitute(Δbasis[j], d))
            end
            b[i] = Float64(substitute(eq, d))
            i += 1
        catch e
            # println("exclude ", x₀)
        end
    end

    l = find_independent_subset(A; abstol)
    A, b, basis = A[l,l], b[l], basis[l]
    q = nice_parameters(A \ b)
    sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x))
end

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x, k=1)
    c = sum(candidates(eq, x))
    D = Differential(x)
    Δc = expand_derivatives(D(c))
    μ = sum(x^i for i=0:k-1)
    # p = expand((1+x) * (c + Δc))
    p = expand(μ * (c + Δc))
    kers = kernelize(p, x)
    kers = expand(kers)
    return [one(x); candidates(kers, x)]
end

function find_independent_subset(A; abstol=1e-5)
    _, R = qr(A)
    abs.(diag(R)) .> abstol
end

function nice_parameters(p; abstol=1e-3)
    c = lcm(collect(1:10)...)
    n = length(p)
    q = Array{Any}(undef, n)
    for i = 1:n
        den = 1
        while den < 10
            if abs(round(p[i]*den) - p[i]*den) < abstol
                a = round(Int, p[i]*den) // den
                q[i] = (denominator(a) == 1 ? numerator(a) : a)
                den = 10
            else
                q[i] = p[i]
            end
            den += 1
        end
    end
    q
end

########################## Transformation Rules ###############################

trig_rule1 = @rule tan(~x) => sin(~x) / cos(~x)
trig_rule2 = @rule sec(~x) => one(~x) / cos(~x)
trig_rule3 = @rule csc(~x) => one(~x) / sin(~x)
trig_rule4 = @rule cot(~x) => cos(~x) / sin(~x)

trig_rules = [trig_rule1, trig_rule2, trig_rule3, trig_rule4]

hyper_rule1 = @rule tanh(~x) => sinh(~x) / cosh(~x)
hyper_rule2 = @rule sech(~x) => one(~x) / cosh(~x)
hyper_rule3 = @rule csch(~x) => one(~x) / sinh(~x)
hyper_rule4 = @rule coth(~x) => cosh(~x) / sinh(~x)

hyper_rules = [hyper_rule1, hyper_rule2, hyper_rule3, hyper_rule4]

misc_rule1 = @rule sqrt(~x) => ^(~x, 0.5)
misc_rule2 = @rule log(~x) => log(abs(~x))
# misc_rule2 = @rule log(~x) => log(~x) * li(~x)

misc_rules = [misc_rule1]

int_rules = [trig_rules; hyper_rules; misc_rules]

apply_integration_rules(eq) = Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq))

########################## Test If a Polynomial? #############################

is_number(x::T) where T<:Integer = true
is_number(x::T) where T<:Real = true
is_number(x::T) where T<:Complex = true
is_number(x::T) where T<:Rational = true
is_number(x) = false

function is_poly(eq, x)
    p = collect_powers(value(eq), x)
    all(is_number, values(p))
end

@parameters a b c d e f

##############################################################################

function read_maxima_test(name)
    fn = open(name, "r")

    for line in readlines(fn)
        println(line)
        if length(line)>1 && line[1] == '['
            line = replace(line, "%pi" => "π")
            line = replace(line, "%e" => "ℯ")

            i = findfirst(isequal(']'), line)
            # try
                l = eval(Meta.parse(line[1:i]))
                p = l[1]
                q = l[4]
                dict = Dict(v => rand(-10:10) for v in get_variables(p) if !isequal(v,x))
                println(dict)
                p = substitute(p, dict)
                q = substitute(q, dict)
                printstyled(p, '\t'; color=:blue)
                printstyled(q, '\t'; color=:green)
                printstyled(integrate(p), '\n'; color=:red)
            # catch e
            # end
        end
    end

    close(fn)
end

@syms x

basic_integrals = [
# Basic Forms
    1,
    x^2,
    4x^3,
# Integrals of Rational Functions
    1/x,
    1/(2x + 5),
    1/(x + 1)^2,
    (x + 3)^3,
    x*(x - 2)^4,
    1/(1 + x^2),
    1/(9 + x^2),
    x/(4 + x^2),
    x^2/(16 + x^2),
    x^3/(1 + x^2),
    1/(x^2 - 5x + 6),
    1/(x^2 + x + 1),
    x/(x + 4)^2,
    x/(x^2 + x + 1),
# Integrals with Roots
    sqrt(x - 2),
    1 / sqrt(x - 1),
    1 / sqrt(x + 1),
    1 / sqrt(4 - x),
    x * sqrt(x - 3),
    sqrt(2x + 5),
    (3x - 1)^1.5,
    x / sqrt(x - 1),
    x / sqrt(x + 1),
    sqrt(x / (4 - x)),
    sqrt(x / (4 + x)),
    x * sqrt(2x + 3),
    sqrt(x*(x+2)),
    sqrt(x^3*(x+3)),
    sqrt(x^2 + 4),
    sqrt(x^2 - 4),
    sqrt(4 - x^2),
    x * sqrt(x^2 + 9),
    x * sqrt(x^2 - 9),
    1 / sqrt(x^2 + 4),
    1 / sqrt(x^2 - 4),
    1 / sqrt(4 - x^2),
    x / sqrt(x^2 + 4),
    x / sqrt(x^2 - 4),
    x / sqrt(4 - x^2),
    x^2 / sqrt(x^2 + 4),
    x^2 / sqrt(x^2 - 4),
    sqrt(x^2 - 5x + 6),
    x * sqrt(x^2 - 5x + 6),
    1 / sqrt(x^2 - 5x + 6),
    1 / (4 + x^2)^1.5,
# Integrals with Logarithms
    log(x),
    x * log(x),
    x^2 * log(x),
    log(2x) / x,
    log(x) / x^2,
    log(2x + 1),
    log(x^2 + 4),
    log(x^2 - 4),
    log(x^2 - 5x + 6),
    x * log(x + 2),
    x * log(9 - 4x^2),
    log(x)^2,
    log(x)^3,
    x * log(x)^2,
    x^2 * log(x)^2,
# Integrals with Exponentials
    exp(x),
    sqrt(x) * exp(x),
    x * exp(x),
    x * exp(3x),
    x^2 * exp(x),
    x^2 * exp(5x),
    x^3 * exp(x),
    x^3 * exp(2x),
    exp(x^2),
    x * exp(x^2),
# Integrals with Trigonometric Functions
    sin(4x),
    sin(x)^2,
    sin(x)^3,
    cos(3x),
    cos(x)^2,
    cos(2x)^3,
    sin(x) * cos(x),
    sin(3x) * cos(5x),
    sin(x)^2 * cos(x),
    sin(3x)^2 * cos(x),
    sin(x) * cos(x)^2,
    sin(x) * cos(5x)^2,
    sin(x)^2 * cos(x),
    sin(x)^2 * cos(x)^2,
    sin(4x)^2 * cos(4x^2),
    tan(x),
    tan(7x),
    tan(x)^2,
    tan(x)^3,
    sec(x),
    sec(x)^2,
    sec(x)^3,
    sec(x) * tan(x),
    sec(x)^2 * tan(x),
    csc(x),
    csc(x)^2,
    csc(x)^3,
    sec(x) * csc(x),
# Products of Trigonometric Functions and Monomials
    x * cos(x),
    x * cos(3x),
    x^2 * cos(x),
    x^2 * cos(5x),
    x * sin(x),
    x * sin(3x),
    x^2 * sin(x),
    x^2 * sin(5x),
    x * cos(x)^2,
    x * sin(x)^2,
    x * tan(x)^2,
    x * sec(x)^2,
# Products of Trigonometric Functions and Exponentials
    exp(x) * sin(x),
    exp(3x) * sin(2x),
    exp(x) * cos(x),
    exp(2x) * cos(7x),
    x * exp(x) * sin(x),
    x * exp(x) * cos(x),
# Integrals of Hyperbolic Functions
    cosh(x),
    exp(x) * cosh(x),
    sinh(3x),
    exp(2x) * sinh(3x),
    tanh(x),
    exp(x) * tanh(x),
    cos(x) * cosh(x),
    cos(x) * sinh(x),
    sin(x) * cosh(x),
    sin(x) * sinh(x),
    sinh(x) * cosh(x),
    sinh(3x) * cosh(5x),
]

function test_integrals()
    for eq in basic_integrals
        printstyled(eq, " =>\t"; color=:green)
        printstyled(integrate(eq), '\n'; color=:red)
    end
end

##################### Special Functions ######################################

function li(x; n=10)
    z = log(abs(x))
    s = sum(z^k / (factorial(k) * k) for k = 1:n)
    return SpecialFunctions.γ + log(z) + s
end
