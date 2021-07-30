using LinearAlgebra

isdependent(eq, x) = !isequal(expand_derivatives(Differential(x)(eq)), 0)

kernelize(eq, x) = x
kernelize(eq::Num, x) = kernelize(value(eq), x)
kernelize(eq::T, x) where T<:Number = one(x)
kernelize(eq::SymbolicUtils.Term, x) = (isdependent(eq,x) ? eq : one(x))
kernelize(eq::SymbolicUtils.Pow, x) = (isdependent(eq,x) && arguments(eq)[2] < 0 ? eq*log(arguments(eq)[1]) : eq)
kernelize(eq::SymbolicUtils.Mul, x) = prod(kernelize(t,x) for t in arguments(eq); init=one(x))
kernelize(eq::SymbolicUtils.Add, x) = sum(kernelize(t,x) for t in arguments(eq) if isdependent(t,x); init=zero(x))

function kernelize(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) || arguments(eq)[2] >= 0
        return eq
    else
        f = factor(arguments(eq)[1])
        return eq * sum(log(f[i]) for i = 1:length(f) if isdependent(f[i],x); init=zero(x))
    end
end

candidates(eq, x) = [eq]
candidates(eq::Num, x) = candidates(value(eq), x)
candidates(eq::SymbolicUtils.Add, x) = unique(∪([candidates(t,x) for t in arguments(eq)]...))

function candidates(eq::SymbolicUtils.Mul, x)
    l = []
    terms = arguments(eq)
    n = length(terms)
    mask = sum(1<<(j-1) for j=1:n if !isdependent(terms[j],x); init=0)

    for i = 1:2^n-1
        if i & mask == 0
            y = prod(terms[j] for j=1:n if ((1<<(j-1)) & i) != 0)
            push!(l, y)
        end
    end
    l
end

function integrate(eq; abstol=1e-5, num_trials=5, num_zeros=5, lo=-5.0, hi=5.0)
    x = var(eq)
    D = Differential(x)

    for i = 1:num_trials
        y = try_integrate(eq; abstol, lo, hi)
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
    end
    @warn "no solution is found"
    nothing
end

function try_integrate(eq; abstol=1e-5, lo=-5.0, hi=5.0)
    eq = apply_integration_rules(eq)
    x = var(eq)

    basis = generate_basis(eq, x)
    # println(basis)
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
    sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0)
end

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x)
    c = sum(candidates(eq, x))
    D = Differential(x)
    Δc = expand_derivatives(D(c))
    p = expand((1+x) * (c + Δc))
    kers = kernelize(p, x)
    kers = expand(kers)    
    return candidates(kers, x)
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

int_rules = [trig_rules; hyper_rules]

apply_integration_rules(eq) = Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq))

##############################################################################

# function quadgk_multi(f, t)
#     n = length(t)-1
#     y = zeros(n)
#     for i = 1:n
#         y[i] = quadgk(f, t[1], t[i+1])[1]
#     end
#     return y
# end
#
# function collect_param_coef(eqs, p)
#     n = length(p)
#     m = length(eqs)
#     A = zeros(m,n)
#     b = zeros(m)
#
#     for (k,eq) in enumerate(eqs)
#         println(eq)
#         d = Dict(p[j] => 0 for j = 1:n)
#         b[k] = Float64(substitute(eq, d))
#         println(b[k])
#
#         for i = 1:n
#             d = Dict(p[j] => (i==j ? 1 : 0) for j = 1:n)
#             A[k,i] = Float64(substitute(eq, d)) - b[k]
#         end
#         println(A[k,:])
#     end
#     return A, -b
# end
