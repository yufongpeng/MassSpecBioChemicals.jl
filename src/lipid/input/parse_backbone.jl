function parse_singlechirality(x)
    if isnothing(x) || isempty(x)
        RSChirality
    elseif x == "(R)-"
        RChirality
    elseif x == "(S)-"
        SChirality
    elseif x == "L-"
        LForm
    elseif x == "D-"
        DForm
    else
        @warn "Invalid chirality"
        RSChirality
    end
end

# specific for glycerol-phospho-aa
function parse_glycerol_rs(pre)
    if isnothing(pre)
        Glycerol{RSChirality}()
    elseif pre == "(R)-"
        Glycerol{RChirality}()
    elseif pre == "(S)-"
        Glycerol{SChirality}()
    else
        @warn "Invalid chirality"
        Glycerol{RSChirality}()
    end
end

function parse_glycerol2_rs(pre)
    isnothing(pre) && return (Glycerol{RSChirality}(), Glycerol{RSChirality}())
    rss = split(replace(pre, "(" => "", ")" => "", " " => "", "-" => ""), ",")
    gly1, gly2 = (Glycerol{RSChirality}(), Glycerol{RSChirality}())
    for rs in rss 
        if rs == "R"
            gly1 = Glycerol{RChirality}()
        elseif rs == "S"
            gly1 = Glycerol{SChirality}()
        elseif rs == "R'"
            gly2 = Glycerol{RChirality}()
        elseif rs == "S'"
            gly2 = Glycerol{SChirality}()
        else
            @warn "Invalid chirality"
        end
    end
    (gly2, gly1)
end

function parse_glycerol3_rs(pre)
    isnothing(pre) && return (Glycerol{RSChirality}(), Glycerol{RSChirality}(), Glycerol{RSChirality}())
    rss = split(replace(pre, "(" => "", ")" => "", " " => "", "-" => ""), ",")
    gly1, gly2, gly3 = (Glycerol{RSChirality}(), Glycerol{RSChirality}(), Glycerol{RSChirality}())
    for rs in rss 
        if rs == "R"
            gly1 = Glycerol{RChirality}()
        elseif rs == "S"
            gly1 = Glycerol{SChirality}()
        elseif rs == "R'"
            gly2 = Glycerol{RChirality}()
        elseif rs == "S'"
            gly2 = Glycerol{SChirality}()
        elseif rs == "R''"
            gly3 = Glycerol{RChirality}()
        elseif rs == "S''"
            gly3 = Glycerol{SChirality}()
        else
            @warn "Invalid chirality"
        end
    end
    (gly3, gly2, gly1)
end

function parse_aa_glycerol_rs(pre, aa = "Ser")
    isnothing(pre) && return (parse_aa(aa), Glycerol{RSChirality}())
    rs, dl = match(r"(?:\(([RS\,']+)\)-)*([DL]-)*", pre)
    if isnothing(dl)
        aaa = parse_aa(aa)
    else
        aaa = parse_aa(string(dl, aa))
    end
    if isnothing(rs)
        gly = Glycerol{RSChirality}()
    else
        rss = split(replace(rs, " " => ""), ",")
        gly = Glycerol{RSChirality}()
        for rs in rss 
            if rs == "R"
                gly = Glycerol{RChirality}()
            elseif rs == "S"
                gly = Glycerol{SChirality}()
            elseif rs == "R'"
                aaa = parse_aa(string("(R)-", aa))
            elseif rs == "S'"
                aaa = parse_aa(string("(S)-", aa))
            else
                @warn "Invalid chirality"
            end
        end
    end
    (aaa, gly)
end