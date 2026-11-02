function Base.getindex(db::LipidDataBase, i)
    LipidDataBaseObject(db, i)
end

function Base.getindex(dbo::LipidDataBaseObject, i)
    dbo.database.records[i][dbo.id]
end

function Base.getindex(db::LipidDataBase, i, j)
    db.database.records[j][i]
end

Base.propertynames(db::LipidDataBase, private = false) =
    private ? (:records, (Symbol(k) for k in keys(db.records))...) : Tuple(Symbol(k) for k in keys(db.records))
Base.propertynames(dbo::LipidDataBaseObject, private = false) =
    private ? (:database, (Symbol(k) for k in keys(dbo.database.records))...) : Tuple(Symbol(k) for k in keys(dbo.database.records))

function Base.getproperty(db::LipidDataBase, i)
    if i == :records
        getfield(db, i)
    else
        getfield(db, :records)[String(i)]
    end
end

function Base.getproperty(dbo::LipidDataBaseObject, i)
    if i == :database || i == :id
        getfield(dbo, i)
    else
        getfield(dbo, :database).records[String(i)][dbo.id]
    end
end

Base.iterate(db::LipidDataBase, state = 1) =
    state > length(db) ? nothing : (LipidDataBaseObject(db, state), state + 1)
Base.length(db::LipidDataBase) = length(first(values(getfield(db, :records))))

Base.keys(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.values(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.pairs(dbo::LipidDataBaseObject) = (k => getproperty(dbo, k) for k in propertynames(dbo))
function Base.iterate(dbo::LipidDataBaseObject, state = 1)
    p = propertynames(dbo)
    state > length(p) ? nothing : (getproperty(dbo, p[state]), state + 1)
end

function Base.show(io::IO, level::AnnotationLevel)
    print(io, @match level begin
        &specieslevel                   => "specieslevel"
        &molecularspecieslevel          => "molecularspecieslevel"
        &phosphatepositionlevel         => "phosphatepositionlevel"
        &snpositionlevel                => "snpositionlevel"
        &dbpositionpartiallevel         => "dbpositionpartiallevel"
        &dbpositionlevel                => "dbpositionlevel"
        &dbconfigpartiallevel           => "dbconfigpartiallevel"
        &dbconfiglevel                  => "dbconfiglevel"
        &structuredefinedpartiallevel   => "structuredefinedpartiallevel"
        &structuredefinedlevel          => "structuredefinedlevel"
        &structurepositionpartiallevel  => "structurepositionpartiallevel"
        &structurepositionlevel         => "structurepositionlevel"
        &structureconfigpartiallevel    => "structureconfigpartiallevel"
        &structureconfiglevel           => "structureconfiglevel"
        &fullstructurelevel             => "fullstructurelevel"
        &completestructurelevel         => "completestructurelevel"
    end)
end

function Base.show(io::IO, level::PassLevel)
    print(io, @match level begin
        &passphosphatepositionlevel     => "passphosphatepositionlevel"
        &passsnpositionlevel            => "passsnpositionlevel"
    end)
end