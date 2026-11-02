config_chain(::Type{<: Lipid}) = UInt8[]
config_chain(::Type{<: Monoradylglycerol}) = [0x02]
config_chain(::Type{<: Diradylglycerol}) = [0x02]
config_chain(::Type{<: Triradylglycerol}) = [0x02]
config_chain(::Type{<: Omodifiedradylglycerol}) = [0x02]
config_chain(::Type{<: Radylglycerophosphate}) = [0x02] 
config_chain(::Type{<: Bisradylglycerophosphate}) = [0x02, 0x04] 
config_chain(::Type{<: Bisradylglycerophosphoglycerol}) = [0x01, 0x03] 
config_chain(::Type{<: SphingoBone}) = [0x02] 
config_chain(::Type{<: MixSphingoBone}) = [0x02] 

check_configuration!(bone, chain::CarbonChain{<: Tuple}) = chain
function check_configuration!(bone::B, chain::CarbonChain{A}; headposition = nothing) where {B, A <: AbstractSPB} 
    free_db = false
    check_db = true
    if ispropertyempty(chain.doublebond) 
        ispropertyempty(chain.substituent) && return chain # cb:0 
        ispropertynumber(chain.substituent) && any(x -> first(x) isa UnknownGroup, chain.substituent) && return chain # cb:0;O\d 
        check_db = false
    elseif ispropertynumber(chain.doublebond) 
        free_db = true 
        check_db = false
    elseif ispropertypartialposition(chain.doublebond)
        free_db = true 
    end
    free_fg = 0 
    free_df = 0
    free_oxo = false
    check_fg  = true
    if ispropertyempty(chain.substituent)
        check_fg = false 
    elseif ispropertynumber(chain.substituent) && any(x -> first(x) isa UnknownGroup, chain.substituent)
        check_fg = false
    elseif ispropertynumber(chain.substituent)
        check_fg = false
        if length(chain.substituent) == 1 && first(first(chain.substituent)) == Oxo()
            free_oxo = true
        else 
            if all(x -> first(x) isa (FunctionalGroup{M, Didehydrogen} where M), chain.substituent) # bridge 
                # check last db
            end
            free_fg = Int(sum(last, chain.substituent))
            free_df = length(filter(x -> first(x) != Oxo() && first(x) != Hydrogen(), chain.substituent))
        end
    elseif ispropertypartialposition(chain.substituent)
        if all(last(x) == Oxo() for x in chain.substituent if first(x) == 0)
            free_oxo = true
        else
            free_fg = count(x -> first(x) == 0, chain.substituent) 
            free_df = length(unique(filter(x -> first(x) == 0 && last(x) != Oxo() && last(x) != Hydrogen(), chain.substituent)))
        end
    end
    if check_db
        db_all = Dict{Int, Any}(Int(k) => v for (k, v) in chain.doublebond)
    else
        db_all = Dict{Int, Any}()
    end
    if check_fg && isnothing(headposition)
        fg_all = Dict{Int, Vector{<: AbstractFunctionalGroup}}()
        for (k, v) in chain.substituent
            f = get!(fg_all, Int(k), AbstractFunctionalGroup[])
            push!(f, v)
        end
    elseif !isnothing(headposition)
        fg_all = Dict{Int, Vector{<: AbstractFunctionalGroup}}()
        for (k, v) in headposition
            f = get!(fg_all, Int(k), AbstractFunctionalGroup[])
            push!(f, v)
        end
    else
        fg_all = Dict{Int, Any}()
    end
    if isnothing(chain.chirality)
        ch_all = Dict{Int, Any}()
    else
        ch_all = Dict{Int, RSSystem}()
        for (k, v) in chain.chirality 
            f = get!(ch_all, Int(k), v)
            f == v && continue 
            if f == RSChirality() 
                ch_all[k] = v
            elseif v == RSChirality() 
                continue 
            else
                throw(ArgumentError("Conflict of chirality at position `$k`."))
            end
        end
        # delete unknown chiral @ db
        for i in keys(db_all)
            if get(ch_all, i, nothing) == RSChirality()
                delete!(ch_all, i)
            end
        end
    end 
    sym = false
    config_check_begin!(bone, chain, fg_all, ch_all, db_all)
    # Assume no symmetry
    for i in 3:Int(chain.carbon) - 1
        config_check_inner!(bone, chain, fg_all, ch_all, db_all, i, sym; free_fg, free_oxo, free_db) 
    end
    sym = config_check_end!(bone, chain, fg_all, ch_all, db_all, sym; free_fg, free_df)
    empty!(chain.chirality)
    for (k, v) in ch_all 
        push!(chain.chirality, UInt8(k) => v)
    end
    sort!(chain.chirality; by = first)
    if check_db 
        empty!(chain.doublebond) 
        for (k, v) in db_all 
            push!(chain.doublebond, UInt8(k) => v)
        end
        sort!(chain.doublebond; by = first)
    end
    chain
end

function check_configuration!(bone::B, chain::CarbonChain{A}) where {B, A <: Radyl} 
    sym = isodd(chain.carbon)
    free_db = false
    check_db = true
    if ispropertyempty(chain.doublebond) 
        ispropertyempty(chain.substituent) && return chain # cb:0;#O, cb:0
        ispropertynumber(chain.substituent) && return chain # cb:0;#FG
        # cb:0;pFG, check sub
        check_db = false
    elseif ispropertynumber(chain.doublebond) 
        sym = false 
        check_db = false
        free_db = true 
    elseif ispropertypartialposition(chain.doublebond)
        sym = false 
        free_db = true 
    end
    free_fg = 0 
    free_df = 0
    free_oxo = false
    check_fg  = true
    if ispropertyempty(chain.substituent) 
        check_fg = false 
    elseif ispropertynumber(chain.substituent) && any(x -> first(x) isa UnknownGroup, chain.substituent)
        check_fg = false
    elseif ispropertynumber(chain.substituent)
        if length(chain.substituent) == 1 && first(first(chain.substituent)) == Oxo()
            free_oxo = true
        else 
            if all(x -> first(x) isa (FunctionalGroup{M, Didehydrogen} where M), chain.substituent) # bridge 
                # check last db
            end
            free_fg = Int(sum(last, chain.substituent))
            free_df = length(filter(x -> first(x) != Oxo() && first(x) != Hydrogen(), chain.substituent))
        end
        check_fg = false
        sym = false
    elseif ispropertypartialposition(chain.substituent)
        sym = false 
        if all(last(x) == Oxo() for x in chain.substituent if first(x) == 0)
            free_oxo = true
        else
            free_fg = Int(sum(last(x) for x in chain.substituent if first(x) == 0)) 
            free_df = length(unique(filter(x -> first(x) == 0 && last(x) != Oxo() && last(x) != Hydrogen(), chain.substituent)))
        end
    end
    if check_db
        db_all = Dict{Int, GeometricConfiguration}(Int(k) => v for (k, v) in chain.doublebond)
    else
        db_all = Dict{Int, Any}()
    end
    if check_fg 
        fg_all = Dict{Int, Vector{<: AbstractFunctionalGroup}}()
        for (k, v) in chain.substituent
            f = get!(fg_all, Int(k), AbstractFunctionalGroup[])
            push!(f, v)
        end
    else
        fg_all = Dict{Int, Any}()
    end
    if isnothing(chain.chirality)
        ch_all = Dict{Int, Any}()
    else
        ch_all = Dict{Int, RSSystem}()
        for (k, v) in chain.chirality 
            f = get!(ch_all, Int(k), v)
            f == v && continue 
            if f == RSChirality() 
                ch_all[k] = v
            elseif v == RSChirality() 
                continue 
            else
                throw(ArgumentError("Conflict of chirality at position `$k`."))
            end
        end
        # delete unknown chiral @ db
        for i in keys(db_all)
            if get(ch_all, i, nothing) == RSChirality()
                delete!(ch_all, i)
            end
        end
    end
    # sym = sym && get(ch_all, Int((chain.carbon + 1) / 2), AChirality()) == RSChirality()
    config_check_begin!(bone, chain, fg_all, ch_all, db_all)
    for i in 2:Int(chain.carbon) - 1
        sym = config_check_inner!(bone, chain, fg_all, ch_all, db_all, i, sym; free_fg, free_oxo, free_db) 
    end
    if sym 
        set_achiral_nowarn!(ch_all, Int((chain.carbon + 1) / 2))
    end
    sym = config_check_end!(bone, chain, fg_all, ch_all, db_all, sym; free_fg, free_df)
    empty!(chain.chirality)
    for (k, v) in ch_all 
        push!(chain.chirality, UInt8(k) => v)
    end
    sort!(chain.chirality; by = first)
    if check_db 
        empty!(chain.doublebond) 
        for (k, v) in db_all 
            push!(chain.doublebond, UInt8(k) => v)
        end
        sort!(chain.doublebond; by = first)
    end
    chain
end

function config_check_begin!(bone::B, chain::CarbonChain{Acyl}, fg_all, ch_all, db_all) where B 
    haskey(db_all, 1) && throw(ArgumentError("Position `1` cannot contain double bond."))
    haskey(fg_all, 1) && throw(ArgumentError("Too many functional groups at position `1`.")) 
    haskey(ch_all, 1) && throw(ArgumentError("Position `1` cannot be a chiral center.")) # warn ?
    # only check fg
end

function config_check_begin!(bone::B, chain::CarbonChain{Alkenyl{C}}, fg_all, ch_all, db_all) where {B, C} 
    haskey(db_all, 1) && throw(ArgumentError("Position `1` cannot contain additional double bond."))
    length(get(fg_all, 1, [])) > 1 && throw(ArgumentError("Too many functional groups at position `1`.")) 
    length(get(fg_all, 2, [])) > 1 && throw(ArgumentError("Too many functional groups at position `2`.")) 
    # check no bridge
    haskey(ch_all, 1) && throw(ArgumentError("Position `1` cannot be a chiral center.")) # warn ?
    haskey(ch_all, 2) && throw(ArgumentError("Position `2` cannot be a chiral center.")) # warn ?
end

function config_check_begin!(bone::B, chain::CarbonChain{Alkyl}, fg_all, ch_all, db_all) where B
    # if haskey(db_all, 1) 
    #     throw(ArgumentError("Position `1` cannot contain double bond for alkyl ether lipid."))
    # end
    if haskey(ch_all, 1) 
        fg_check2(fg_all, 1)
        v = vcat(v, fg_all[1])
        ch_check3!(ch_all, 1, v)
    end
    if haskey(db_all, 1) 
        dbw = (db_all[1] == EConfiguration() || db_all[1] == ZConfiguration())
        f = filter!(!=(Hydrogen()), get(fg_all, 1, [dehydrogengroup(bone)]))
        if isempty(f)
            dbw && @warn "Change position `i` to no E/Z configuration."
            db_all[1] = NoEZConfiguration()
        elseif length(f) == 2 && ischemicalequal(f...)
            dbw && @warn "Change position `i` to no E/Z configuration."
            db_all[1] = NoEZConfiguration()
        end
        if haskey(fg_all, 2) && chain.carbon >= 3
            # check fg 
            f = first(fg_all[2])
            # Overwritten explicit EZ
            if ishydrocarbon(f) && ncarbon(f) == chain.carbon - 2 && all(x -> !haskey(fg_all, x) && !haskey(db_all, x), 3:chain.carbon) 
                dbw && @warn "Change position `1` to no E/Z configuration."
                db_all[1] = NoEZConfiguration()
            end
        end
    end
end 

function config_check_begin!(bone::B, chain::CarbonChain{<: AbstractSPB}, fg_all, ch_all, db_all) where B
    if haskey(db_all, 1) 
        throw(ArgumentError("Position `1` cannot contain double bond."))
    end
    if haskey(db_all, 2) 
        throw(ArgumentError("Position `2` cannot contain double bond."))
    end
    if haskey(ch_all, 2) 
        first(fg_all[2]) == Hydrogen() || throw(ArgumentError("Position `2` cannot contain funtional group."))
    end
    if haskey(ch_all, 1) 
        fg_check3(fg_all, 1)
        ch_check3!(ch_all, 1, fg_all[1])
    end
end 

function config_check_end!(bone::B, chain::CarbonChain{Acyl}, fg_all, ch_all, db_all, sym; free_fg = 0, free_df = 0) where B
    v = (isnothing(bone) || bone == Dihydrogen()) ? Oxo() : XLinkedFunctionalGroup(CarboxylicLinkage(), dehydrogengroup(bone))
    _config_check_end!(bone, chain, fg_all, ch_all, db_all; free_fg, free_df)
    if sym 
        if bone == HydrogenOxide()
            sym = in(CarboxylicAcidGroup(), get(fg_all, Int(chain.carbon), []))
        else
            lfg = get(fg_all, Int(chain.carbon), [])
            filter!(!=(Hydrogen()), lfg)
            sym = any(x -> v == x, lfg)
        end
        f = get(ch_all, 1, AChirality())
        g = get(ch_all, Int(chain.carbon), AChirality())
        sym = sym && f == g && f != RSChirality()
    end
    sym
end

function config_check_end!(bone::B, chain::CarbonChain{Alkenyl{C}}, fg_all, ch_all, db_all, sym; free_fg = 0, free_df = 0) where {B, C} 
    # only check fg
    _config_check_end!(bone, chain, fg_all, ch_all, db_all; free_fg, free_df)
    if sym 
        v = filter(!=(Hydrogen()), vcat(isnothing(bone) ? [] : [dehydrogengroup(bone)], get(fg_all, 1, [])))
        lfg = filter(!=(Hydrogen()), get(fg_all, Int(chain.carbon), []))
        sym = get(db_all, Int(chain.carbon) - 1, nothing) == C() && all(x -> in(x, lfg), v)
        f = get(ch_all, 1, AChirality())
        g = get(ch_all, Int(chain.carbon), AChirality())
        sym = sym && f == g && f != RSChirality()
    end
    sym
end

function config_check_end!(bone::B, chain::CarbonChain{Alkyl}, fg_all, ch_all, db_all, sym; free_fg = 0, free_df = 0) where B
    v = isnothing(bone) ? AbstractFunctionalGroup[] : AbstractFunctionalGroup[dehydrogengroup(bone)]
    _config_check_end!(bone, chain, fg_all, ch_all, db_all; free_fg, free_df)
    if sym 
        # check last 
        v = filter(!=(Hydrogen()), v)
        lfg = filter(!=(Hydrogen()), get(fg_all, Int(chain.carbon), []))
        lfg = vcat(lfg, filter(x -> x == Epoxy() || x == Peroxy(), get(fg_all, Int(chain.carbon) - 1, [])))
        sym = !haskey(db_all, Int(chain.carbon) - 1) && all(x -> in(x, lfg), v)
        f = get(ch_all, 1, AChirality())
        g = get(ch_all, Int(chain.carbon), AChirality())
        sym = sym && f == g && f != RSChirality()
    end
    sym
end

function config_check_end!(bone::B, chain::CarbonChain{<: AbstractSPB}, fg_all, ch_all, db_all, sym; free_fg = 0, free_df = 0) where B
    _config_check_end!(bone, chain, fg_all, ch_all, db_all; free_fg, free_df)
    sym
end

function _config_check_end!(bone::B, chain::CarbonChain, fg_all, ch_all, db_all; free_fg = 0, free_df = 0) where B 
    i = Int(chain.carbon)
    if haskey(ch_all, i) 
        v = fg_all[i]
        fg_check3(fg_all, i)
        ch_check3!(ch_all, i, v; free_fg, free_df)
        if ch_all[i] != AChirality() && haskey(db_all, i - 1) 
            throw(ArgumentError("Position `$(i - 1)` cannot contain double bond."))
        end
    end
end

function fg_check2(fg_all, i)
    sum(nlinkage, fg_all[i]) > 2 && throw(ArgumentError("Too many functional groups at position `$i`.")) 
    # length(v) > 2 && throw(ArgumentError("Too many functional groups at position `$i`.")) 
    # if length(v) > 1
    #     n = length(findall(==(Oxo()), v))
    #     n > 0 && throw(ArgumentError("Too many functional groups at position `$i`.")) 
    # end
end

function fg_check3(fg_all, i)
    sum(nlinkage, fg_all[i]) > 3 && throw(ArgumentError("Too many functional groups at position `$i`.")) 
    # v = fg_all[i]
    # println(v)
    # length(v) > 3 && throw(ArgumentError("Too many functional groups at position `$i`.")) 
    # if length(v) > 1
    #     n = length(findall(==(Oxo()), v))
    #     if n > 1
    #         throw(ArgumentError("Multiple oxo cannot locate at the same position `$i`."))
    #     elseif n == 1 && length(v) > 2
    #         throw(ArgumentError("Too many functional groups at position `$i`.")) 
    #     end
    # end
end

function ch_check3!(ch_all, i, v; free_fg = 0, free_df = 0)
    n = count(==(Didehydrogen()), v)
    n > 1 && return set_achiral_nowarn!(ch_all, i)
    n == 1 && (free_fg < 2 || free_df < 2) && return set_achiral_nowarn!(ch_all, i)
    if length(v) == 1
        # free_df == 1 && check == free
        if free_fg == 0
            return set_achiral_nowarn!(ch_all, i)
        end
    elseif length(v) == 2 
        if ischemicalequal(v[2], v[1])
            return set_achiral_nowarn!(ch_all, i)
        end
    else
        (
            ischemicalequal(v[1], v[2]) ||
            ischemicalequal(v[2], v[3]) ||
            ischemicalequal(v[1], v[3])
        ) && return set_achiral_nowarn!(ch_all, i)

    end
end

function set_achiral!(ch_all, i)
    v = get!(ch_all, i, AChirality())
    ch_all[i] = AChirality()
    v == AChirality() || @warn "Change Position `$i` to achiral."
end

function set_achiral_warn!(ch_all, i)
    ch_all[i] = AChirality()
    @warn "Change Position `$i` to achiral."
end

function set_achiral_nowarn!(ch_all, i)
    ch_all[i] = AChirality()
end

function config_check_inner!(bone::B, chain::CarbonChain{A}, fg_all, ch_all, db_all, i, sym; free_fg = 0, free_oxo = false, free_db = false) where {B, A} 
    hc = B <: Dihydrogen && A <: Alkyl
    check_cf = free_fg == 0 && !free_oxo && !free_db
    haskey(fg_all, i) && fg_check2(fg_all, i) # check fg oxo => only one fg 
    if haskey(ch_all, i) && ch_all[i] != RSChirality() && ch_all[i] != AChirality()
        # check db
        if haskey(db_all, i) || haskey(db_all, i - 1)
            throw(ArgumentError("Position `$(i - 1)` and `$i` cannot contain double bond."))
        end
        c = ch_all[i]
        v = get(fg_all, i, AbstractFunctionalGroup[])
        if length(v) == 2 && ischemicalequal(v...)
            set_achiral_nowarn!(ch_all, i)
        elseif length(v) == 1 && free_fg == 0 && first(v) == Hydrogen() && all(x -> !in(x, [Epoxy(), Peroxy()]), get(fg_all, i - 1, []))
            set_achiral_nowarn!(ch_all, i)
        elseif c == RSChirality()
            for f in v 
                if check_cf && ishydrocarbon(f) 
                    # not allow modification and db on fg
                    if ncarbon(f) == i - 1 && hc && all(x -> !haskey(fg_all, x) && !haskey(db_all, x), 1:i - 1) 
                        set_achiral_nowarn!(ch_all, i)
                    elseif ncarbon(f) == chain.carbon - i && all(x -> !haskey(fg_all, x) && !haskey(db_all, x), i + 1:chain.carbon) 
                        set_achiral_nowarn!(ch_all, i)
                    end
                end
            end                        
        end
    end
    if haskey(db_all, i) 
        dbw = (db_all[i] == EConfiguration() || db_all[i] == ZConfiguration())
        if haskey(fg_all, i)
            # check fg 
            length(fg_all[i]) > 1 && throw(ArgumentError("Position `$i` can not have two or more substitutions."))
            haskey(ch_all, i) && (ch_all[i] = AChirality())
            f = first(fg_all[i])
            # Overwritten explicit EZ
            if check_cf && ishydrocarbon(f) 
                # not allow modification and db on fg
                if ncarbon(f) == i - 1 && hc && all(x -> !haskey(fg_all, x) && !haskey(db_all, x), 1:i - 1) 
                    dbw && @warn "Change position `$i` to no E/Z configuration."
                    db_all[i] = NoEZConfiguration()
                end
            end
        end
        if i + 1 == chain.carbon
            if haskey(fg_all, i + 1) 
                length(fg_all[i + 1]) > 2 && throw(ArgumentError("Position `$(i + 1)` can not have three or more substitutions."))
                haskey(ch_all, i + 1) && (ch_all[i + 1] = AChirality())
                # Overwritten explicit EZ
                if length(fg_all[i + 1]) == 2 && ischemicalequal(fg_all[i + 1]...)
                    dbw && @warn "Change position `$i` to no E/Z configuration."
                    db_all[i] = NoEZConfiguration()
                end
            elseif free_fg == 0 
                dbw && @warn "Change position `$i` to no E/Z configuration."
                db_all[i] = NoEZConfiguration()
            end
        elseif haskey(fg_all, i + 1) 
            length(fg_all[i + 1]) > 1 && throw(ArgumentError("Position `$(i + 1)` can not have two or more substitutions."))
            haskey(ch_all, i + 1) && (ch_all[i + 1] = AChirality())
            f = first(fg_all[i + 1])
            # Overwritten explicit EZ
            if check_cf && ishydrocarbon(f) 
                # not allow modification and db on fg
                if ncarbon(f) == chain.carbon - i - 1 && all(x -> !haskey(fg_all, x) && !haskey(db_all, x), i + 2:chain.carbon) 
                    dbw && @warn "Change position `$i` to no E/Z configuration."
                    db_all[i] = NoEZConfiguration()
                end
            end
        end
    end
    if sym && i > (chain.carbon + 1) ÷ 2
        # check other side
        v = filter(!=(Hydrogen()), get(fg_all, Int(chain.carbon) - i + 1, []))
        lfg = filter(!=(Hydrogen()), get(fg_all, i, []))
        lfg = vcat(lfg, filter(x -> x isa Epoxy() || x isa Peroxy(), get(fg_all, i - i, [])))
        sym = all(x -> in(x, lfg), v)
        f = get(ch_all, i, AChirality())
        g = get(ch_all, Int(chain.carbon) - i + 1, AChirality())
        sym = sym && f == g && f != RSChirality()
        f = get(db_all, i - 1, NoEZConfiguration())
        g = get(db_all, Int(chain.carbon) - i + 1, NoEZConfiguration())
        sym = sym && f == g && f != EZConfiguration()
        # println(i => sym)
    end
    sym
end