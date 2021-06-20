module SymPoly

using Symbolics, SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, solve_for, derivative
using SymbolicUtils.Rewriters
using SymbolicUtils.Code

include("poly.jl")

export Poly, extract_coef, gen_poly, degree, terms, var, leading, cont, prim
export poly, polyas, sym, update_order

include("arith.jl")

export divide_poly, derivative, integrate, rationalize, to_monic, from_monic, integer_poly

include("roots.jl")

export solve_newton, find_roots

include("factors.jl")

export FactoredPoly, factors, factor_schubert_kronecker, decompose, recompose, factor

include("tests.jl")

end # module
