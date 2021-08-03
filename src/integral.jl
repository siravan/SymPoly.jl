using LinearAlgebra
# using SpecialFunctions

"""
    isdependent returns true if eq is dependent on x
"""
isdependent(eq, x) = !isequal(expand_derivatives(Differential(x)(eq)), 0)

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x)
    Δeq = expand_derivatives(Differential(x)(eq))
    kers = expand(eq + Δeq)
    return [one(x); candidates(kers, x)]
end

"""
    candidates returns a list of candidate expressions to form the integration
    basis
"""
candidates(eq, x) = isdependent(eq,x) ? [eq] : []

candidates(eq::Num, x) = candidates(value(eq), x)

# the candidates of an Add is the union of the candidates of the terms
# ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
candidates(eq::SymbolicUtils.Add, x) = unique(∪([candidates(t,x) for t in arguments(eq)]...))

# the candidates of a Mul is the outer product of the candidates of the terms
# d(uv)/dx = u dv/dx + v + du/dx
function candidates(eq::SymbolicUtils.Mul, x)
    terms = [candidates(q,x) for q in arguments(eq)]
    n = length(terms)

    l = Any[one(x)]

    for j = 1:n
        m = length(l)
        for t in terms[j]
            for k = 1:m
                push!(l, l[k]*t)
            end
        end
    end

    unique(l[2:end])    # removing the initial 1
end

# the candidates of a Pow encode different integration rules
function candidates(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) return [one(x)] end

    p = arguments(eq)[1]    # eq = p ^ k
    k = arguments(eq)[2]

    if k < 0 && k ≈ round(k)
        return candidate_pow_minus(p, k, x)
    elseif k ≈ 0.5 || k ≈ -0.5
        # ∫ √f(x) dx = ... + c * log(df/dx + √f) if deg(f) == 2
        Δ = expand_derivatives(Differential(x)(p))
        return [[p^k, p^(k+1)]; log(abs(0.5*Δ + sqrt(p)))]
    end

    # ∫ p^k dp = c * p^(k+1)
    return [p^k, p^(k+1)]
end

function candidate_pow_minus(p, k, x)
    q = to_poly(p, x)   # q is a DynamicPolynomials Polynomial
                        # the reason is to be able to use find_roots
    if q == nothing return [p^k, p^(k+1), log(p)] end

    q = 0 + 1*q         # take care of monomials and constants
    r, s = find_roots(q, var(q))
    s = s[1:2:end]
    r = nice_parameters(r)
    s = Complex.(nice_parameters(real.(s)), nice_parameters(imag(s)))

    # ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
    # q = sum(log(x - u) for u in r; init=zero(x)) +
    #     sum(atan((x - real(u))/imag(u)) for u in s; init=zero(x)) +
    #     sum(log(x^2 - 2*real(u)*x + abs2(u)) for u in s; init=zero(x))

    q₁ = [[log(x - u) for u in r];
          [atan((x - real(u))/imag(u)) for u in s];
          [log(x^2 - 2*real(u)*x + abs2(u)) for u in s]
         ]

    q₂ = [[(x - u)^-1 for u in r];
          [(x^2 - 2*real(u)*x + abs2(u))^-1 for u in s];
          [x*(x^2 - 2*real(u)*x + abs2(u))^-1 for u in s]
         ]

    # return [[p^k, p^(k+1)]; candidates(q₁, x)]
    if k ≈ -1
        return [[p^k]; q₁; q₂]
    else
        return [[p^k, p^(k+1)]; q₁; q₂]
    end
end

###############################################################################

"""
    integrate is the main entry point

    input:
    ------
    eq: a Symbolics expression to integrate
    abstol: the desired tolerance
    num_steps: the number of different steps with expanding basis to be tried
    num_trials: the number of trials in each step (no changes to the basis)
    lo and hi: the range used to generate random values of x (the independent variable)
    show_basis: if true, the basis is printed

    output:
    -------
    solved, unsolved

    a pair of expressions, solved is the solved integral and unsolved is the residual unsolved
    portion of the input
"""
function integrate(eq; abstol=1e-6, num_steps=3, num_trials=2, lo=-5.0, hi=5.0, show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false)
    x = var(eq)
    eq = expand(eq)

    # eq is a constant
    if x == nothing
        @syms 𝑥
        return 𝑥 * eq, 0
    end

    # check if eq is a rational function
    # if so, we perform a partial-fraction decomposition first (the first part of the Hermite's method)

    # q = to_rational(eq, x)
    # if q != nothing
    #     eq = q
    # end

    s₁, u₁ = integrate_sum(eq, x; abstol, num_trials, num_steps, lo, hi, show_basis, opt, bypass)

    if isequal(u₁, 0)
        return s₁, u₁
    else
        s₂, u₂ = try_integration_by_parts(u₁, x; abstol, num_trials, num_steps, lo, hi, show_basis, opt)
        return s₁ + s₂, u₂
    end
end

"""
    ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
"""
function integrate_sum(eq::SymbolicUtils.Add, x; abstol=1e-6, num_steps=3, num_trials=2, lo=-5.0, hi=5.0, show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false)
    if bypass
        return integrate_term(eq, x; abstol, num_steps, num_trials, lo, hi, show_basis, opt)
    else
        solved = 0
        unsolved = 0

        for p in arguments(eq)
            s, u = integrate_term(p, x; abstol, num_steps, num_trials, lo, hi, show_basis, opt)
            solved += s
            unsolved += u
        end

        return solved, unsolved
    end
end

function integrate_sum(eq, x; abstol=1e-6, num_steps=3, num_trials=2, lo=-5.0, hi=5.0, show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false)
    integrate_term(eq, x; abstol, num_steps, num_trials, lo, hi, show_basis, opt)
end

function integrate_term(eq, x; abstol=1e-6, num_steps=3, num_trials=2, lo=-5.0, hi=5.0, show_basis=false, opt = STLSQ(exp.(-10:1:0)))
    eq₁ = apply_integration_rules(eq)
    basis = generate_basis(eq₁, x)

    D = Differential(x)

    for i = 1:num_steps
        basis = unique([basis; basis*x])
        Δbasis = [expand_derivatives(D(f)) for f in basis]
        if show_basis println(basis) end

        for j = 1:num_trials
            y, ϵ = try_integrate(eq, x, basis, Δbasis; abstol, lo, hi, opt)
            if ϵ < abstol return y, 0 end
        end
    end
    # @warn "no solution is found"
    0, eq
end

rms(x) = sqrt(sum(x.^2) / length(x))

"""
    the core of the randomized parameter-fitting algorithm

    `try_integrate` tries to find a linear combination of the basis, whose
    derivative is equal to eq

    output
    -------
    integral, error
"""
function try_integrate(eq, x, basis, Δbasis; abstol=1e-6, lo=-5.0, hi=5.0, attemp_ratio=5, opt = STLSQ(exp.(-10:1:0)))
    n = length(basis)
    # A is an nxn matrix holding the values of the fragments at n random points
    # b hold the value of the input function at those points
    A = zeros(n, n)
    b = zeros(n)

    i = 1
    k = 1

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
        if k > attemp_ratio*n return nothing, 1e6 end
        k += 1
    end

    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, b, basis = A[l,l], b[l], basis[l]

    if det(A) ≈ 0 return nothing, 1e6 end

    # q₀ = A \ b
    q₀ = Optimize.init(opt, A, b)
    @views Optimize.sparse_regression!(q₀, A, permutedims(b)', opt, maxiter = 1000)
    q = nice_parameters(q₀)
    ϵ = rms(A * q - b)
    sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x)), ϵ
end

"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol=1e-5)
    _, R = qr(A)
    abs.(diag(R)) .> abstol
end

"""
    converts float to int or small rational numbers
"""
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

function try_integration_by_parts(eq, x; abstol=1e-6, num_steps=3, num_trials=2, lo=-5.0, hi=5.0, show_basis=false, opt = STLSQ(exp.(-10:1:0)))
    f = factors(eq, x)
    if length(f) <= 2 return zero(x), eq end

    D = Differential(x)

    for u in f
        v′ = eq / u
        if !is_number(u) && !is_number(v′)
            v, r = integrate_term(v′, x; abstol, num_steps, num_trials, lo, hi, show_basis, opt)
            if isequal(r, 0)
                uv = expand_derivatives(v*D(u))
                s, r = integrate_sum(uv, x; abstol, num_steps, num_trials, lo, hi, show_basis, opt)
                if isequal(r, 0)
                    return expand(u*v - s), 0
                end
            end
        end
    end

    return zero(x), eq
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

misc_rules = [misc_rule1]

int_rules = [trig_rules; hyper_rules; misc_rules]

apply_integration_rules(eq) = Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq))

########################## Convert to a Polynomial? #############################

@polyvar 𝑦

is_number(x::T) where T<:Integer = true
is_number(x::T) where T<:Real = true
is_number(x::T) where T<:Complex = true
is_number(x::T) where T<:Rational = true
is_number(x) = false

function is_poly(eq, x)
    p = collect_powers(value(eq), x)
    all(is_number, values(p))
end

function to_poly(p::SymbolicUtils.Add, x)
    l = [to_poly(t,x) for t in arguments(p)]
    if any(x->x==nothing, l) return nothing end
    sum(l; init=0)
end

function to_poly(p::SymbolicUtils.Mul, x)
    l = [to_poly(t,x) for t in arguments(p)]
    if any(x->x==nothing, l) return nothing end
    prod(l; init=1)
end

function to_poly(p::SymbolicUtils.Pow, x)
    if !isequal(arguments(p)[1], x) return nothing end
    return 𝑦^arguments(p)[2]
end

to_poly(p::SymbolicUtils.Sym, x) = isequal(p, x) ? 𝑦 : nothing
to_poly(p::SymbolicUtils.Term, x) = nothing
to_poly(p, x) = is_number(p) ? p : nothing

function to_rational(p::SymbolicUtils.Mul, x)
    P = one(x)
    Q = one(x)

    for t in arguments(p)
        y = to_poly(t, x)
        if y != nothing
            P *= y
        else
            y = to_poly(1/t, x)
            if y != nothing
                Q *= y
            else
                return nothing
            end
        end
    end

    if is_number(P) || is_number(Q) return nothing end
    return sym(factor_rational(P, Q), x)
end

to_rational(p, x) = nothing

is_multiple_x(p::SymbolicUtils.Add, x) = all(is_multiple_x(t,x) for t in arguments(p))
is_multiple_x(p::SymbolicUtils.Mul, x) = any(is_multiple_x(t,x) for t in arguments(p))
is_multiple_x(p::SymbolicUtils.Pow, x) = isequal(arguments(p)[1], x) && arguments(p)[2] >= 1
is_multiple_x(p::SymbolicUtils.Sym, x) = isequal(p, x)
is_multiple_x(p, x) = false

factors(eq, x) = isdependent(eq, x) ? [one(x), eq] : [one(x)]

function factors(eq::SymbolicUtils.Pow, x)
    p, k = arguments(eq)
    [p^(i*sign(k)) for i=0:abs(k)]
end

function factors(eq::SymbolicUtils.Mul, x)
    terms = [factors(q,x) for q in arguments(eq)]
    n = length(terms)

    l = Any[one(x)]

    for j = 1:n
        m = length(l)
        for t in terms[j]
            for k = 1:m
                push!(l, l[k]*t)
            end
        end
    end

    unique(l)
end

##############################################################################

@syms x β

"""
    a list of basic standard integral tests
    based on http://integral-table.com/ with modifications
"""
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
    x^3 * sin(x),
    x^4 * cos(2x),
    sin(x)^2 * cos(x)^3,
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
# Misc
    exp(x) / (1 + exp(x)),
    cos(exp(x)) * sin(exp(x)) * exp(x),
    cos(exp(x))^2 * sin(exp(x)) * exp(x),
    1 / (x*log(x)),
    (log(x - 1) + (x - 1)^-1) * log(x),
    1 / (exp(x) - 1),
    1 / (exp(x) + 5),
    sqrt(x)*log(x),
    log(log(x)) / x,
    x^3 * exp(x^2),
    cos(log(x)),
    x * cos(x) * exp(x),
    log(x - 1)^2,
    1 / (exp(2x) - 1),
    exp(x) / (exp(2x) - 1),
    x / (exp(2x) - 1),
# derivative-divide examples (Lamangna 7.10.2)
    exp(x) * exp(exp(x)),
    exp(sqrt(x)) / sqrt(x),
    log(log(x)) / (x*log(x)),
    log(cos(x)) * tan(x),
# rothstein-Trager examples (Lamangna 7.10.9)
    1 / (x^3 - x),
    1 / (x^3 + 1),
    1 / (x^2 - 8),
    (x + 1) / (x^2 + 1),
    x / (x^4 - 4),
    x^3 / (x^4 + 1),
    1 / (x^4 + 1),
# bypass = true
    β,      # turn of bypass = true
    exp(x)/x - exp(x)/x^2,
    cos(x)/x - sin(x)/x^2,
    1/log(x) - 1/log(x)^2,
]

function test_integrals()
    bypass = false

    for eq in basic_integrals
        if isequal(eq, β)
            printstyled("**** bypass on ****\n"; color=:red)
            bypass = true
        else
            printstyled(eq, " =>\t"; color=:green)
            solved, unsolved = integrate(eq; bypass)
            printstyled(solved; color=:white)
            if isequal(unsolved, 0)
                println()
            else
                printstyled(" + ∫ ", unsolved, '\n'; color=:red)
            end
        end
    end
end

##################### Special Functions ######################################

# """
#     logarithmic integral
# """
# function li(x; n=10)
#     z = log(abs(x))
#     s = sum(z^k / (factorial(k) * k) for k = 1:n)
#     return SpecialFunctions.γ + log(z) + s
# end
