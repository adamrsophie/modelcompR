# Dumps the coefficient tables of the stabilized Runge-Kutta methods of
# OrdinaryDiffEqStabilizedRK (ROCK2, ROCK4, ESERK4, ESERK5, SERK2) and of
# the stabilized SDE methods of StochasticDiffEqROCK (SROCK2, SKSROCK,
# TangXiaoSROCK2, KomBurSROCK2, SROCKC2) as exact hexadecimal floats, one
# field per line: cache, field, dimensions, values. tools/stabilized-gen.R
# turns the dump into R/sysdata.rda.
#
# Run from the repository root:
#   julia --project=tools/julia tools/stabilized-dump.jl tools/stabilized-tables.tsv
#   Rscript tools/stabilized-gen.R tools/stabilized-tables.tsv
using OrdinaryDiffEqStabilizedRK, StochasticDiffEq, Printf
const SR = OrdinaryDiffEqStabilizedRK
const RK = Base.require(Base.PkgId(Base.UUID("db241ea8-0e6b-4abc-8f2d-1adff2294fd9"), "StochasticDiffEqROCK"))
hex(x) = @sprintf("%a", Float64(x))
function dump(io, name, c)
    for fn in fieldnames(typeof(c))
        v = getfield(c, fn)
        fn in (:zprev,) && continue
        if v isa AbstractVector{<:Tuple}
            # a vector of coefficient tuples: one row per entry
            v = permutedims(reduce(hcat, [collect(x) for x in v]))
        end
        if v isa Union{AbstractArray{<:Real}, Tuple{Vararg{Real}}}
            a = v isa Tuple ? collect(v) : v
            dims = join(size(a), "x")
            println(io, name, "\t", fn, "\t", dims, "\t", join(hex.(vec(a)), ","))
        elseif v isa Real
            println(io, name, "\t", fn, "\t1\t", hex(v))
        end
    end
end
open(ARGS[1], "w") do io
    T = Float64
    dump(io, "ROCK2", SR.ROCK2ConstantCache(T, T, zeros(1)))
    dump(io, "ROCK4", SR.ROCK4ConstantCache(T, T, zeros(1)))
    dump(io, "ESERK4", SR.ESERK4ConstantCache(T, zeros(1)))
    dump(io, "ESERK5", SR.ESERK5ConstantCache(T, zeros(1)))
    dump(io, "SERK2", SR.SERK2ConstantCache(T, zeros(1)))
    for nm in (:SROCK2, :SKSROCK, :TangXiaoSROCK2, :KomBurSROCK2, :SROCKC2, :SROCK1, :SROCKEM)
        c = getfield(RK, Symbol(nm, :ConstantCache)){T}(zeros(1))
        dump(io, string(nm), c)
    end
end
