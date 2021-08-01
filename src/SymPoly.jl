module SymPoly

using Primes

using MultivariatePolynomials
using DynamicPolynomials

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters
# using SymbolicUtils.Code

using Symbolics

include("symbolics.jl")

export poly, sym

include("factoredpoly.jl")
include("wrapper.jl")

export wrap, unwrap, var, leading, cont, prim
export deg, derivative
export rationalize, unrationalize
export gcd_extended

include("roots.jl")

export solve_newton, find_roots

include("factorization.jl")

export FactoredPoly, factors, factor_schubert_kronecker
export decompose, factor, expand_frac, power, frac

include("modular.jl")
include("roundabout.jl")

include("roots_comb.jl")

include("integral.jl")

export integrate, generate_basis, test_integrals

#############################################################################

function Primes.factor(p; method=:roots_SSP)
    if method == :schubert_kronecker
        return factor_schubert_kronecker(p)
    # elseif method == :roundabout
    #    return factor_roundabout(p)
    elseif method == :roots_comb
        return factor_roots_comb(p)
    elseif method == :roots_SSP
        return factor_roots_SSP(p)
    else
        error("undefined factorization method")
    end
end

Primes.factor(p, q) = factor_rational(p, q)
Primes.factor(r::RationalPoly) = factor_rational(r)

end # module
