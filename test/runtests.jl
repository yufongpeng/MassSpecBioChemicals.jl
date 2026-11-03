using MassSpecBioChemicals
using MassSpecBioChemicals.Lipids
const MBL = MassSpecBioChemicals.Lipids
using Test, JSON3

function test_annotationlevel(x)
    tal = MBL.annotationlevel(x.object; partial = true, additional = true, pass = true)
    length(tal) == length(x.annotationlevel) && all(t -> in(t, x.annotationlevel), tal)
end

test_show(x) = show(IOBuffer(), x)
macro test_noerror(x)
    x = esc(x)
    return quote
        try 
            $x
            true
        catch e
            false
        end
    end
end

@testset "MassSpecBioChemicals.jl" begin
    @testset "MassSpecBioChemicals.Lipids" begin
        # use json
        test_lipid_js = JSON3.read(joinpath("data", "test_lipid.json"))
        @testset "stdin/stdout" begin
            global test_lipid = Dict{UnionAll, Dict}()
            for (c, s) in test_lipid_js
                @testset "$c" begin
                    dict = Dict{String, Any}()
                    for (l, a) in s
                        l = string(l)
                        lipid = parse_chemical(LipidParser(), l)
                        @test chemicalname(lipid) == first(a)
                        push!(dict, l => (object = lipid, annotationlevel = [eval(Meta.parse(aa)) for aa in last(a)]))
                    end
                    push!(test_lipid, eval(c) => dict)
                end
            end
        end
        @testset "base" begin 
            for (c, d) in test_lipid
                @testset "$c" begin
                    for (k, x) in d
                        @test copy(x.object) == x.object
                        @test hash(x.object) == hash(deepcopy(x.object))
                    end
                end
            end
        end
        @testset "annotationlevel" begin
            for (c, d) in test_lipid
                @testset "$c" begin
                    for (k, x) in d
                        @test test_annotationlevel(x)
                    end
                end
            end
        end
        @testset "SMILES" begin 
            for (c, d) in test_lipid
                @testset "$c" begin
                    for (k, v) in d 
                        @test @test_noerror chemicalsmiles(v.object; onlycarbonchain = true)
                    end
                end
            end
        end
    end
end
