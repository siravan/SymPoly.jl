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
                push!(l, numerator(p) / (denominator(p))^k)
            end
        end
        return length(l)==0 ? 0 : sum(l)
    else
        return prod(map(v -> first(v)^last(v), factors(f)); init=1)
    end
end

function show_cont(io, c, i)
    if i == 1
        if c < 0
            print(io, "-")
        end
    else
        print(io, c > 0 ? " + " : " - ")
    end
end

function Base.show(io::IO, f::FactoredPoly)
    i = 1
    if f.rational
        for w in factors(f)
            p, k = first(w), last(w)

            if p isa AbstractPolynomial
                c = cont(p)
                show_cont(io, c, i)
                if isone(c)
                    print(io, prim(p))
                else
                    print(io, abs(c), " * (", prim(p), ")")
                end
            elseif p isa RationalPoly
                c = cont(numerator(p)) // cont(denominator(p))
                show_cont(io, c, i)
                c = abs(c)
                if isone(c)
                    print(io, '(', prim(numerator(p)), ") / (", prim(denominator(p)), ")")
                else
                    print(io, c, " * (", prim(numerator(p)), ") / (", prim(denominator(p)), ")")
                end

                if k != 1
                    print(io, "^", k)
                end
            else
                print(io, p)
            end
            i += 1
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
        if v isa AbstractPolynomialLike
            add_factor!(h, sym(v, r), k)
        elseif v isa RationalPoly
            cn, n = integer_poly(numerator(v))
            cd, d = integer_poly(denominator(v))
            c = cn ÷ cd
            add_factor!(h, sym(c*n, r) / sym(d, r)^k, 1)
        else
            add_factor!(h, v, k)
        end
    end
    h
end

SymbolicUtils.simplify(f::FactoredPoly) = prod(first(w)^last(w) for w in factors(f); init=1)

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
