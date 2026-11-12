# Exported names of the SciML sublibraries modelcompR ports, read from their
# source (no package is loaded), for tools/ledger.R.
#
# Run from the repository root:
#   julia tools/upstream_exports_src.jl
# Reads the monorepo sources of OrdinaryDiffEq (lib/ holds the solver
# sublibraries) and Optimization (lib/ holds the solver wrappers) from the
# Julia depot and writes tools/upstream_exports_src.tsv (package, name).
# A `@reexport using M` is listed as the line `package  @reexport M`.
function latest(pkg)
    d = joinpath(first(DEPOT_PATH), "packages", pkg)
    vs = filter(isdir, joinpath.(d, readdir(d)))
    vs[argmax(mtime.(vs))]
end
function exports_of(dir)
    out = String[]
    for (root, _, files) in walkdir(joinpath(dir, "src")), f in files
        endswith(f, ".jl") || continue
        ex = try
            Meta.parseall(read(joinpath(root, f), String))
        catch
            continue
        end
        walk(e) = if e isa Expr
            if e.head in (:export, :public)
                append!(out, string.(filter(x -> x isa Symbol, e.args)))
            elseif e.head == :macrocall && length(e.args) >= 3 && e.args[1] == Symbol("@reexport")
                u = e.args[end]
                if u isa Expr && u.head == :using
                    for a in u.args
                        a isa Expr && a.head == :. && push!(out, "@reexport " * join(string.(a.args), "."))
                        a isa Expr && a.head == :(:) && push!(out, "@reexport " * join(string.(a.args[1].args), "."))
                    end
                end
            else
                foreach(walk, e.args)
            end
        end
        walk(ex)
    end
    unique(out)
end
open(joinpath(@__DIR__, "upstream_exports_src.tsv"), "w") do io
    println(io, "package\tname")
    for mono in ("OrdinaryDiffEq", "Optimization")
        lib = joinpath(latest(mono), "lib")
        for sub in sort(readdir(lib))
            isdir(joinpath(lib, sub, "src")) || continue
            for n in sort(exports_of(joinpath(lib, sub)))
                println(io, sub, "\t", n)
            end
        end
    end
    for pkg in ("SciMLOperators", "DiffEqNoiseProcess", "OptimizationBase")
        for n in sort(exports_of(latest(pkg)))
            println(io, pkg, "\t", n)
        end
    end
end
