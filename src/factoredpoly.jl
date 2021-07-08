using Primes

struct FactoredPoly
    rational::Bool
    factors::Array{Pair{Any,Int}}

    FactoredPoly(rational::Bool=false) = new(rational, Array{Pair{Any,Int}}[])
    FactoredPoly(x) = new(false, [x => 1])
end

factors(f::FactoredPoly) = f.factors

# function add_factor!(f::FactoredPoly, w::Wrapper, power::Int=1)
#     push!(f.factors, w.p => power)
# end

function add_factor!(f::FactoredPoly, p, power::Int=1)
    push!(f.factors, p => power)
end

function combine!(f::FactoredPoly, g::FactoredPoly)
    f.rational != g.rational && error("cannot combine incompatible FactoredPolys")
    append!(f.factors, g.factors)
end

Base.length(f::FactoredPoly) = length(f.factors)
Base.getindex(f::FactoredPoly, k::Integer) = first(f.factors[k])
Base.getindex(f::FactoredPoly, r::UnitRange) = first.(f.factors[r])
power(f::FactoredPoly, k::Integer) = last(f.factors[k])

function poly(f::FactoredPoly)
    if f.rational
        l = []
        for w in factors(f)
            p, k = first(w), last(w)
            if p isa AbstractPolynomial
                push!(l, p)
            else
                push!(l, numerator(p) / rationalize(denominator(p))^k)
            end
        end
        return length(l)==0 ? 0 : sum(l)
    else
        return prod(map(v -> first(v)^last(v), factors(f)); init=1)
    end
end

function Base.show(io::IO, f::FactoredPoly)
    i = 1
    if f.rational
        for w in factors(f)
            i > 1 && print(io, " + ")
            i += 1
            p, k = first(w), last(w)

            if p isa AbstractPolynomial
                c, q = integer_poly(p)
                if isone(c)
                    print(io, q)
                else
                    print(io, c, " * (", q, ")")
                end
            elseif p isa RationalPoly
                cn, n = integer_poly(numerator(p))
                cd, d = integer_poly(denominator(p))
                c = cn ÷ cd

                if !isone(c)
                    print(io, c, " * ")
                end

                print(io, '(', n, ") / (", d, ")")

                if k != 1
                    print(io, "^", k)
                end
            else
                print(io, p)
            end
        end
    else
        for w in factors(f)
            p, k = first(w), last(w)
            i > 1 && print(io, " * ")
            i += 1
            if k == 1
                print(io, '(', p, ")")
            else
                print(io, '(', p, ")^", k)
            end
        end
    end
end

function sym(f::FactoredPoly, r)
    h = FactoredPoly(f.rational)
    for w in factors(f)
        v, k = first(w), last(w)
        if v isa AbstractPolynomial
            add_factor!(h, sym(v, r), k)
        elseif v isa RationalPoly
            cn, n = integer_poly(numerator(v))
            cd, d = integer_poly(denominator(v))
            c = cn ÷ cd
            add_factor!(h, sym(c*n, r) / sym(d, r)^k, 1)
        end
    end
    h
end

"""
    convert rational coefficients with a denominator of 1 to integer
"""
function unrationalize(f::FactoredPoly)
    h = FactoredPoly()
    for w in factors(f)
        p, k = first(w), last(w)
        add_factor!(h, unrationalize(p), k)
    end
    h
end
