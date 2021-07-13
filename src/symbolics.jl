@syms 𝑥
@polyvar 𝑦

function poly(eq, v::Pair)
    α, β = value(first(v)), last(v)
    terms = collect_powers(eq, α)

    if !all(isone, denominator.(collect(keys(terms))))
        error("fractional power not supported")
    end

    sum([val*β^numerator(k) for (k,val) in terms]; init=0)
end

function poly(eq, β)
    xs = get_variables(eq)

    if length(xs) > 1
        error("more than one implicit variable.")
    elseif length(xs) == 0
        return nothing
    else
        return poly(eq, xs[1] => β)
    end
end

poly(eq) = poly(eq, 𝑦)

function sym(p::AbstractPolynomialLike, v::Pair)
    β, α = first(v), value(last(v))
    sum([α^maxdegree(t,β)*c for (t,c) in zip(terms(p), coefficients(p))]; init=0)
end

sym(p::AbstractPolynomialLike, α) = sym(p, var(p) => α)
sym(p::AbstractPolynomialLike) = sym(p, var(p) => 𝑥)

###############################################################################

# pox (power-of-x) is a symbolic function to keep track of the powers of x
# pox(k,n) means k*x^n
@syms pox(k, n)

is_pox(x) = istree(x) && operation(x)==pox
is_not_pox(x) = !is_pox(x)

get_coef(p) = is_pox(p) ? arguments(p)[1] : p
get_power(p) = is_pox(p) ? arguments(p)[2] : 0

replace_x(eq, x) = substitute(eq, Dict(x => pox(1,1)))

iscomplex(x) = x isa Complex

count_rule1 = @rule ^(pox(~k, ~n1), ~n2) => isequal(~k,1) ? pox(1, ~n1 * ~n2) : pox(^(~k,~n2), ~n1 * ~n2)
count_rule2 = @rule pox(~k1, ~n1) * pox(~k2, ~n2) => pox(~k1 * ~k2, ~n1 + ~n2)
count_rule3 = @acrule pox(~k, ~n) * ~u::is_not_pox => pox(~k * ~u, ~n)

"""
    collect_powers separates the powers of x in eq (a polynomial) and returns
    a dictionary of power => term
"""
function collect_powers(eq, x)
    eq = expand(expand_derivatives(eq))
    eq = replace_x(eq, x)
    #eq = Prewalk(PassThrough(count_rule1))(eq)
    eq = Fixpoint(Prewalk(PassThrough(Chain([count_rule1, count_rule2, count_rule3]))))(eq)

    if !istree(eq)
        return Dict{Any, Any}(0 => eq)
    elseif is_pox(eq)
        return Dict{Any, Any}(get_power(eq) => get_coef(eq))
    else
        eqs = Dict{Any, Any}()
        for term in arguments(eq)
            n = get_power(term)
            if haskey(eqs, n)
                eqs[n] = eqs[n] + get_coef(term)
            else
                eqs[n] = get_coef(term)
            end
        end

        return eqs
    end
end
