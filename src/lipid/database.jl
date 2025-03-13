# DATABASE[id] -> lazy object
# DATABASE[id][prop] -> lazy getindex prop
# DATABASE[id].prop -> lazy to prop
# DATABASE[id, prop]
# DATABASE.prop
# DATABASE.prop[id]
abstract type LipidDataBase end
abstract type LIPIDMAPS <: LipidDataBase end
struct LMSD <: LIPIDMAPS 
    records::Dict{String, String}
end
struct LipidDataBaseObject{L}
    database::LipidDatabase{L}
    id::Int
end

function init_lipiddatabase!(::Type{LMSD})
    if to_update_LMSD()
        data = take!(Downloads.download("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip", IOBuffer()))
        # write sdf 
        write(joinpath(@__DIR__(), "database", "LMSD.sdf.zip"), data)  
        open(joinpath(@__DIR__(), "database", "status.txt"), "a") do f 
            write(f, "\n", "LMSD\t", Dates.now())
        end
    else
        data = read(joinpath(@__DIR__(), "database", "LMSD.sdf.zip"))
    end
    archive = ZipReader(data)
    sdfs = zip_readentry(archive, "structures.sdf", String)
    # extract data directly from ZipReader buffer
    # split 0x24 0x24 0x24 0x24 0x0a
    sdfs = filter!(!isempty, split(sdfs, "\$\$\$\$\n"))
    mols = Dict(
        "LM_ID" => repeat([""], length(sdfs)), 
        "NAME" => repeat([""], length(sdfs)), 
        "SYSTEMATIC_NAME" => repeat([""], length(sdfs)), 
        "SYNONYMS" => repeat([""], length(sdfs)), 
        "CATEGORY" => repeat([""], length(sdfs)), 
        "MAIN_CLASS" => repeat([""], length(sdfs)), 
        "SUB_CLASS" => repeat([""], length(sdfs)), 
        "CLASS_LEVEL4" => repeat([""], length(sdfs)), 
        "EXACT_MASS" => repeat([""], length(sdfs)), 
        "FORMULA" => repeat([""], length(sdfs)), 
        "INCHI_KEY" => repeat([""], length(sdfs)), 
        "INCHI" => repeat([""], length(sdfs)), 
        "SMILES" => repeat([""], length(sdfs)), 
        "ABBREVIATION" => repeat([""], length(sdfs)), 
        "PUBCHEM_CID" => repeat([""], length(sdfs)), 
        "HMDB_ID" => repeat([""], length(sdfs)), 
        "KEGG_ID" => repeat([""], length(sdfs)), 
        "CHEBI_ID" => repeat([""], length(sdfs)), 
        "LIPIDBANK_ID" => repeat([""], length(sdfs)), 
        "PLANTFA_ID" => repeat([""], length(sdfs)), 
        "SWISSLIPIDS_ID" => repeat([""], length(sdfs)),
        "SDF" => sdfs
    )
    for (i, s) in enumerate(sdfs)
        for (k, v) in eachmatch(r"> <([^>]*)>\n(.*)", s)
            mols[k][i] = v 
        end
    end
    # sdfs = map(read_sdf_prop, sdfs)
    # mols = Dict(k => get.(sdfs, Ref(k), "") for k in [
    #         "LM_ID", 
    #         "NAME", 
    #         "SYSTEMATIC_NAME", 
    #         "SYNONYMS", 
    #         "CATEGORY", 
    #         "MAIN_CLASS", 
    #         "SUB_CLASS", 
    #         "CLASS_LEVEL4", 
    #         "EXACT_MASS", 
    #         "FORMULA", 
    #         "INCHI_KEY", 
    #         "INCHI", 
    #         "SMILES", 
    #         "ABBREVIATION", 
    #         "PUBCHEM_CID", 
    #         "HMDB_ID", 
    #         "KEGG_ID", 
    #         "CHEBI_ID", 
    #         "LIPIDBANK_ID", 
    #         "PLANTFA_ID", 
    #         "SWISSLIPIDS_ID"
    #     ]
    # )
    # push!(mols, "SDF" => sdf)
    # mols = [sdftomol(IOBuffer(x)) for x in filter!(!isempty, split(sdfs, "\$\$\$\$\n"))]
    LMSD(mols)
end

# function read_sdf_prop(m)
#     Dict{String, String}(eachmatch(r"> <([^>]*)>\n(.*)", m))
# end

function to_update_LMSD()
    mkpath(joinpath(@__DIR__(), "database"))
    dir = readdir(joinpath(@__DIR__(), "database"))
    "status.txt" in dir || return true
    st = readlines(joinpath(@__DIR__(), "database", "status.txt"))
    i = findlast(x -> startswith(x, "LMSD\t"), st)
    return isnothing(i) || begin
        Dates.DateTime(
            first(match(
                r"LMSD (\d\d\d\d-\d\d-\d\d)", 
                String(take!(Downloads.download("https://www.lipidmaps.org/databases/lmsd/download", IOBuffer())))
            )), dateformat"yyyy-mm-dd\THH:MM:SS.s"
        ) > Dates.DateTime(last(split(st[i], "\t")), dateformat"yyyy-mm-dd\THH:MM:SS.s")
    end
end

const DATABASE_LMSD = init_lipiddatabase!(LMSD)

const DATABASE_LIPID = Dict{String, LipidDatabase}(
    "LMSD" => DATABASE_LMSD
)

