# module SymPoly

using Reexport

using MultivariatePolynomials
using DynamicPolynomials

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters
# using SymbolicUtils.Code

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

# end # module
