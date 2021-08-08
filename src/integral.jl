using LinearAlgebra
# using SpecialFunctions

"""
    isdependent returns true if eq is dependent on x
"""
isdependent(eq, x) = !isequal(expand_derivatives(Differential(x)(eq)), 0)

Base.signbit(z::Complex{T}) where T<:Number = signbit(real(z))

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x, h=[])
    Δeq = expand_derivatives(Differential(x)(eq))
    kers = expand(eq + Δeq)
    # return [one(x); candidates(kers, x); h]
    return [one(x); candidates(kers, x); h]
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
        # return [[p^k, p^(k+1)]; log(abs(0.5*Δ + sqrt(p)))]
        return [[p^k, p^(k+1)]; log(0.5*Δ + sqrt(p))]
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
    r = nice_parameter.(r)
    s = nice_parameter.(s)
    # s = Complex.(nice_parameters(real.(s)), nice_parameters(imag(s)))

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
function integrate(eq; kwargs...)
    x = var(eq)
    integrate(eq, x; kwargs...)
end

function integrate(eq, x; abstol=1e-6, num_steps=2, num_trials=10, radius=5.0,
                   show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false,
                   attempt_ratio=5, symbolic=false, bypart=true, max_basis=110,
                   verbose=false)
    eq = expand(eq)

    # eq is a constant
    if x == nothing
        @syms 𝑥
        return 𝑥 * eq, 0, 0
    end

    # check if eq is a rational function
    # if so, we perform a partial-fraction decomposition first (the first part of the Hermite's method)

    # q = to_rational(eq, x)
    # if q != nothing
    #     eq = q
    # end

    s₁, u₁, ϵ = integrate_sum(eq, x; bypass, abstol, num_trials, num_steps,
                              radius, show_basis, opt, attempt_ratio, symbolic, max_basis, verbose)

    if isequal(u₁, 0) || !bypart
        return s₁, u₁, ϵ
    else
        s₂, u₂, ϵ = try_integration_by_parts(u₁, x; abstol, num_trials, num_steps,
                                             radius, show_basis, opt, attempt_ratio,
                                             symbolic, max_basis, verbose)
        return s₁ + s₂, u₂, ϵ
    end
end

"""
    ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
"""
function integrate_sum(eq::SymbolicUtils.Add, x; bypass=false, kwargs...)
    if bypass
        return integrate_term(eq, x; kwargs...)
    else
        solved = 0
        unsolved = 0
        ϵ₀ = 0

        for p in arguments(eq)
            s, u, ϵ = integrate_term(p, x; kwargs...)
            solved += s
            unsolved += u
            ϵ₀ = max(ϵ₀, ϵ)
        end

        return solved, unsolved, ϵ₀
    end
end

function integrate_sum(eq, x; kwargs...)
    integrate_term(eq, x; kwargs...)
end

function accept_solution(eq, x, sol; abstol=1e-6)
    try
        Δ = substitute(expand_derivatives(Differential(x)(sol)-eq), Dict(x => Complex(rand())))
        return abs(Δ) < abstol
    catch e
        #
    end
    return false
end

function integrate_term(eq, x; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis, radius =
        args[:abstol], args[:num_steps], args[:num_trials], args[:show_basis],
        args[:symbolic], args[:verbose], args[:max_basis], args[:radius]

    # note that the order of the operations is important!
    # first, collecing hints, then applying transformation rules, and finally finding the basis.
    h = collect_hints(eq, x)
    eq = apply_integration_rules(eq)
    basis = generate_basis(eq, x, h)

    if verbose printstyled("|β| = ", length(basis), ". "; color=:yellow) end
    if length(basis) > max_basis return 0, eq, Inf end

    D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    for i = 1:num_steps
        basis = unique([basis; basis*x])
        Δbasis = [expand_derivatives(D(f)) for f in basis]
        if show_basis println(basis); println(Δbasis) end

        if symbolic
            y, ϵ = try_symbolic(Float64, eq, x, basis, Δbasis; kwargs...)
            if !isequal(y, 0) && accept_solution(eq, x, y; abstol)
                if verbose printstyled("$i, symbolic\n"; color=:yellow) end
                return y, 0, 0
            end
        end

        for j = 1:num_trials
            y, ϵ = try_integrate(Complex{Float64}, eq, x, basis, Δbasis, radius*rand(); kwargs...)
            if ϵ < abstol && accept_solution(eq, x, y; abstol)
                if verbose printstyled("$i, $j\n"; color=:yellow) end
                return y, 0, ϵ
            else
                ϵ₀ = min(ϵ, ϵ₀)
                y₀ = y
            end
        end
    end

    if accept_solution(eq, x, y₀; abstol)
        if verbose printstyled("rescue\n"; color=:yellow) end
        return y₀, 0, ϵ₀
    else
        return 0, eq, ϵ₀
    end
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
function try_integrate(T, eq, x, basis, Δbasis, radius; kwargs...)
    args = Dict(kwargs)
    abstol, opt, attempt_ratio = args[:abstol], args[:opt], args[:attempt_ratio]

    n = length(basis)
    # A is an nxn matrix holding the values of the fragments at n random points
    # b hold the value of the input function at those points
    A = zeros(T, (n, n))
    b = zeros(T, n)

    i = 1
    k = 1

    while i <= n
        if T isa Complex
            x₀ = T(radius*randn()*cis(2π*rand()))
        else
            x₀ = T(radius*rand())
        end
        d = Dict(x => x₀)
        try
            for j = 1:n
                A[i, j] = T(substitute(Δbasis[j], d))
            end
            b[i] = T(substitute(eq, d))
            i += 1
        catch e
            println("basis matrix error: ", e)
        end
        if k > attempt_ratio*n return nothing, 1e6 end
        k += 1
    end

    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, b, basis, Δbasis, n = A[l,l], b[l], basis[l], Δbasis[l], sum(l)

    if det(A) ≈ 0 return nothing, 1e6 end

    coefs = ones(T, n)
    for j = 1:n
        coefs[j] = coef(Δbasis[j], x)
        A[:,j] /= coefs[j]
    end

    # q₀ = A \ b
    q₀ = Optimize.init(opt, A, b)
    @views Optimize.sparse_regression!(q₀, A, permutedims(b)', opt, maxiter = 1000)
    ϵ = rms(A * q₀ - b)
    q = nice_parameter.(q₀ ./ coefs)
    sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x)), abs(ϵ)
end

function find_independent_basis(T, eq, x, Δbasis; kwargs...)
    args = Dict(kwargs)
    abstol, opt, radius, attempt_ratio = args[:abstol], args[:opt], args[:radius], args[:attempt_ratio]

    n = length(Δbasis)
    A = zeros(T, (n, n))

    i = 1
    k = 1

    while i <= n
        if T isa Complex
            x₀ = T(radius*randn()*cis(2π*rand()))
        else
            x₀ = T(radius*rand())
        end
        d = Dict(x => x₀)
        try
            for j = 1:n
                A[i, j] = T(substitute(Δbasis[j], d))
            end
            i += 1
        catch e
            println(e)
        end
        if k > attempt_ratio*n return nothing, 1e6 end
        k += 1
    end

    # find a linearly independent subset of the basis
    find_independent_subset(A; abstol)
end


"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol=1e-3)
    Q, R = qr(A)
    abs.(diag(R)) .> abstol
end

function find_independent_subset2(A; abstol=1e-3)
    n = size(A, 1)
    l = BitVector(undef, n)
    for i = 1:n
        l[i] = 1
        if abs(det(A[l,l])) < abstol
            l[i] = 0
        end
    end
    l
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
                q[i] = Float64(p[i])
            end
            den += 1
        end
    end
    q
end

function nice_parameter(u::T; abstol=1e-3, M=10) where T<:Real
    c = lcm(collect(1:M)...)
    for den = 1:M
        try
            if abs(round(u*den) - u*den) < abstol
                a = round(Int, u*den) // den
                return (denominator(a) == 1 ? numerator(a) : a)
            end
        catch e
        end
    end
    return u
end

function nice_parameter(u::Complex{T}; abstol=1e-3, M=10) where T<:Real
    α = nice_parameter(real(u))
    β = nice_parameter(imag(u))
    return β ≈ 0 ? α : Complex(α, β)
end

function try_integration_by_parts(eq, x; kwargs...)
    f = factors(eq, x)
    if length(f) <= 2 return zero(x), eq, Inf end

    D = Differential(x)
    ϵ₀ = Inf

    for u in f
        v′ = eq / u
        if !is_number(u) && !is_number(v′)
            v, r, ϵ = integrate_term(v′, x; kwargs...)
            if isequal(r, 0)
                uv = expand_derivatives(v*D(u))
                s, r, ϵ = integrate_sum(uv, x; kwargs...)
                if isequal(r, 0)
                    return expand(u*v - s), 0, ϵ
                else
                    zero(x), eq, ϵ
                    # ϵ₀ = min(ϵ, ϵ₀)
                end
            end
        end
    end

    return zero(x), eq, ϵ₀
end

var_index(v) = istree(v) && Symbol(operation(v)) == :getindex ? arguments(v)[2] : -1

function try_symbolic(T, eq, x, basis, Δbasis; kwargs...)
    n = length(basis)
    @syms θ[1:n]
    D = Differential(x)
    Δeq = expand(sum(θ[j]*Δbasis[j] for j=1:n) - eq)

    terms = collect_terms(Δeq, x)
    eqs = collect(values(terms))
    # eqs = filter(p->length(get_variables(p)) > 0, eqs)

    sol = solve_symbolic(eqs)

    for i = 1:n
        if !haskey(sol, θ[i])
            sol[θ[i]] = 0
        end
    end

    p = substitute(expand(sum(θ[j]*basis[j] for j=1:n)), sol)
    p, 0
end

mutable struct Fragment
    eq
    lhs
end

function solve_symbolic(eqs)
    n = length(eqs)
    solved = Set()
    unfinished = true
    frags = Fragment[]
    k = 1

    while unfinished
        unfinished = false
        for eq in eqs
            δf = [v for v in get_variables(eq) if v ∉ solved]
            if length(δf) == 1
                push!(solved, δf[1])
                push!(frags, Fragment(eq, δf[1]))
                unfinished = true
            end
        end
    end

    sol = Dict()
    for f in frags
        eq = substitute(f.eq, sol)
        u = Symbolics.solve_for(eq ~ 0, f.lhs)
        sol[f.lhs] = nice_parameter(u)
    end

    sys = []
    vars = Set()

    for eq in eqs
        δf = [v for v in get_variables(eq) if v ∉ solved]
        if length(δf) > 1
            for v in δf
                push!(vars, v)
            end
            q = substitute(eq, sol)
            push!(sys, q)
        end
    end

    sys = unique(sys)
    vars = [v for v in vars]

    if !isempty(vars) && length(vars) == length(sys)
        try
            vals = nice_parameters.(Symbolics.solve_for(sys .~ 0, vars))
            # vals = Symbolics.solve_for(sys .~ 0, vars)
            for (v,u) in zip(vars, vals)
                sol[v] = u
            end
        catch e
            # println("from symbolic: ", e)
        end
    end

    sol
end

########################## Transformation Rules ###############################

trig_rule1 = @rule tan(~x) => sin(~x) / cos(~x)
trig_rule2 = @rule sec(~x) => one(~x) / cos(~x)
# trig_rule2 = @rule sec(~x) => (tan(~x)/cos(~x) + 1/cos(~x)^2) / (tan(~x) + 1/cos(~x))
trig_rule3 = @rule csc(~x) => one(~x) / sin(~x)
# trig_rule3 = @rule csc(~x) => (cot(~x)/sin(~x) + 1/sin(~x)^2) / (cot(~x) + 1/sin(~x))
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
# int_rules = misc_rules

apply_integration_rules(eq) = Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq))


function U(u...)
    u = map(x -> x isa AbstractArray ? x : [], u)
    return union(u...)
end

hints(eq::SymbolicUtils.Add, x, h) = map(t->hints(t,x,h), arguments(eq))
hints(eq::SymbolicUtils.Mul, x, h) = map(t->hints(t,x,h), arguments(eq))
hints(eq::SymbolicUtils.Pow, x, h) = hints(arguments(eq)[1],x,h)

function hints(eq::SymbolicUtils.Term, x, h)
    s = Symbol(operation(eq))
    u = arguments(eq)[1]

    if s == :sec
        push!(h, log(1/cos(u) + sin(u)/cos(u)))
    elseif s == :csc
        push!(h, log(1/sin(u) - cos(u)/sin(u)))
    elseif s == :tan
        push!(h, log(cos(u)))
    elseif s == :cot
        push!(h, log(sin(u)))
    elseif s == :tanh
        push!(h, log(cosh(u)))
    end
end

function hints(eq, x, h)
end

function collect_hints(eq, x)
    h = []
    hints(eq, x, h)
    h
end

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

coef(eq::SymbolicUtils.Mul, x) = prod(t for t in arguments(eq) if !isdependent(t,x); init=1)
coef(eq::SymbolicUtils.Add, x) = minimum(abs(coef(t,x)) for t in arguments(eq))
coef(eq, x) = 1

function collect_terms(eq::SymbolicUtils.Add, x)
    d = Dict{Any,Any}(1 => 0)
    for t in arguments(eq)
        if isdependent(t, x)
            for s in collect_terms(t, x)
                v, c = first(s), last(s)
                if haskey(d, v)
                    d[v] += c
                else
                    d[v] = c
                end
            end
        else
            d[1] += t
        end
    end
    d
end

function collect_terms(eq, x)
    c = coef(eq, x)
    return Dict(eq/c => c)
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
    sec(x) * tan(x),
    sec(x)^2 * tan(x),
    csc(x),
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
    sin(log(x)),
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

function test_integrals(; symbolic=false, num_trials=3, verbose=false, bypart=true)
    bypass = false
    misses = []
    D = Differential(x)

    for eq in basic_integrals
        if isequal(eq, β)
            printstyled("**** bypass on ****\n"; color=:red)
            bypass = true
        else
            printstyled(eq, " =>\t"; color=:green)
            solved, unsolved = integrate(eq; bypass, symbolic, num_trials, verbose, bypart)
            printstyled(solved; color=:white)
            if isequal(unsolved, 0)
                println()
            else
                printstyled(" + ∫ ", unsolved, '\n'; color=:red)
                push!(misses, eq)
            end
        end
    end

    n = length(misses)
    if n > 0 println("**** missess (n=$n) *****") end
    for eq in misses
        printstyled(eq, '\n'; color=:red)
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
