# DATABASE[id] -> lazy object
# DATABASE[id][prop] -> lazy getindex prop
# DATABASE[id].prop -> lazy to prop
# DATABASE[id, prop]
# DATABASE.prop
# DATABASE.prop[id]
abstract type LipidDataBase end
abstract type LIPIDMAPS <: LipidDataBase end
struct LMSD <: LIPIDMAPS 
    records::Dict{String, Vector{String}}
end
struct LipidDataBaseObject{L <: LipidDataBase}
    database::L
    id::Int
end

function init_lipiddatabase!(::Type{LMSD}; period = Dates.Month(0))
    if to_update_LMSD(period)
        data = take!(Downloads.download("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip", IOBuffer()))
        # write sdf 
        write(joinpath(@__DIR__(), "database", "LMSD.sdf.zip"), data)  
        open(joinpath(@__DIR__(), "database", "status.txt"), "a") do f 
            write(f, "\n", "LMSD\t", string(Dates.now()))
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
    LMSD(mols)
end

function to_update_LMSD(period::DatePeriod)
    mkpath(joinpath(@__DIR__(), "database"))
    dir = readdir(joinpath(@__DIR__(), "database"))
    "status.txt" in dir || return true
    st = readlines(joinpath(@__DIR__(), "database", "status.txt"))
    i = findlast(x -> startswith(x, "LMSD\t"), st)
    isnothing(i) && return true 
    ld = Dates.DateTime(last(split(st[i], "\t")), dateformat"yyyy-mm-dd\THH:MM:SS.s")
    if Dates.now() > ld + period 
        Dates.DateTime(
            first(match(
                r"LMSD (\d\d\d\d-\d\d-\d\d)", 
                String(take!(Downloads.download("https://www.lipidmaps.org/databases/lmsd/download", IOBuffer())))
            )), dateformat"yyyy-mm-dd\THH:MM:SS.s"
        ) > ld
    else 
        false
    end
end

const DATABASE_LMSD = init_lipiddatabase!(LMSD; period = Dates.Month(1))

const DATABASE_LIPID = Dict{String, LipidDataBase}(
    "LMSD" => DATABASE_LMSD
)
