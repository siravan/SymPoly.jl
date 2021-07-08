using Primes

# based on https://discourse.julialang.org/t/arithmetic-modulo-primes/23895/3

struct ℤₚ{p} <: Number
   val::Int64

    function ℤₚ{n}(a) where {n}
        u = mod(a,n)
        new(u > n ÷ 2 ? u - n : u)
    end
end

val(x::ℤₚ{n}) where n = x.val

ℤₚ{n}(x::ℤₚ{n}) where n = x

Base.promote(x::ℤₚ{n}, y::Integer) where {n}=(x,ℤₚ{n}(y))
Base.promote(y::Integer, x::ℤₚ{n}) where {n}=(ℤₚ{n}(y),x)

Base.zero(::Type{ℤₚ{n}}) where {n} = ℤₚ{n}(0)
Base.one(::Type{ℤₚ{n}}) where {n} = ℤₚ{n}(1)

Base.:(==)(x::ℤₚ,y::ℤₚ) = (x.val==y.val)
Base.:(==)(x::ℤₚ{n}, k::Integer) where n = (val(x) == k)
Base.:(==)(k::Integer, x::ℤₚ) = (x == k)

Base.:+(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)+y.val)
Base.:*(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)*y.val)
Base.:-(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(Int(x.val)-y.val)
Base.:-(x::ℤₚ{n}) where {n} = ℤₚ{n}(-Int(x.val))

Base.:/(x::ℤₚ{n}, y::ℤₚ{n}) where n = x * inv(y)
Base.:÷(x::ℤₚ{n}, y::ℤₚ{n}) where n = x * inv(y)
Base.:÷(x::ℤₚ{n}, y::Integer) where n = ℤₚ{n}(Int(x.val)*invmod(y,n))

Base.inv(x::ℤₚ{n}) where {n} = ℤₚ{n}(invmod(x.val,n))
Base.real(x::ℤₚ{n}) where {n} = x.val
Base.abs(x::ℤₚ{n}) where {n} = abs(x.val)
Base.gcd(x::ℤₚ{n}, y::ℤₚ{n}) where {n} = ℤₚ{n}(gcd(x.val, y.val))
Base.gcd(x::ℤₚ{n}, y::Integer) where {n} = gcd(x.val, ℤₚ{n}(y))
Base.gcd(x::Integer, y::ℤₚ{n}) where {n} = gcd(ℤₚ{n}(x), y)

function Base.show(io::IO, m::ℤₚ{n}) where n
    if get(io,:limit, false)
        sub = Dict(zip("0123456789,()","₀₁₂₃₄₅₆₇₈₉‚₍₎"))
        print(io, m.val, map(x->sub[x],repr(n)))
   else
        print(io,"ℤₚ{$n}($(m.val))")
   end
end

###############################################################################

modular(n, v::AbstractArray) = [ℤₚ{n}(x) for x in v]

function modular(n::Integer, p::AbstractPolynomialLike)
     polynomial(modular(n, coefficients(p)), terms(p))
end

# function modular(n::Integer, p::Polynomial{true, T}) where T<:Integer
#     polynomial(modular(n, coefficients(p)), terms(p))
# end
#
# function modular(n::Integer, p::Polynomial{true, ℤₚ{m}}) where m
#     if n == m
#         return p
#     else
#         modular(n, demodular(p))
#     end
# end
#
function demodular(p::AbstractPolynomialLike) where n
    polynomial(val.(coefficients(p)), terms(p))
end

demodular(p::Polynomial{true, T}) where T<:Integer = p
