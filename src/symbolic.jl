# Simplify an expression by applying a few algebraic rules
simplify(e::Number) = e
simplify(e::Symbol) = e

function simplify(e::Expr)
    @assert e.head == :call

    op = e.args[1]
    nargs = length(e.args) - 1
    @assert nargs > 0
    args = e.args[2:end]
    sargs = map(simplify, args)

    
    if op == :* && any(x->x == 0, sargs)
        return 0.0
    end

    if op == :*
        sargs = filter(x->x != 1, sargs)
        if length(sargs) == 1
            return sargs[1]
        elseif length(sargs) == 0
            return 1.0
        end
    end

    if op == :+
        sargs = filter(x->x != 0, sargs)
        if length(sargs) == 1
            return sargs[1]
        elseif length(sargs) == 0
            return 0.0
        end
    end

    if op == :- && nargs == 2 && sargs[2] == 0
        return sargs[1]
    end

    if op == :/
        if sargs[1] == 0 && sargs[1] != 0
            return 0.0
        elseif sargs[2] == 1
            return sargs[1]
        end
    end

    if op == :^
        if sargs[2] == 0
            return 1.0
        elseif sargs[1] == 1
            return 1.0
        elseif sargs[1] == 0
            return 0.0
        elseif sargs[2] == 1
            return sargs[1]
        end
    end
    
    se = expr(:call, { op, sargs... })

    if all(x->isa(x, Number), sargs)
        eval(se)
    else
        se
    end
end

# Computes symbolic derivation of e w.r.t. v
# Should not contain leads or lags on variables
deriv(e::Number, v::Symbol) = 0.0
deriv(e::Symbol, v::Symbol) = e == v ? 1.0 : 0.0
deriv(e::Expr, v::Symbol) = simplify(deriv1(e, v))

function deriv1(e::Expr, v::Symbol)
    @assert e.head == :call

    op = e.args[1]
    
    nargs = length(e.args) - 1
    @assert nargs > 0
    args = e.args[2:end]
    dargs = map(x->deriv(x, v), args)
    if op == :+
        expr(:call, { :+, dargs... })
    elseif op == :*
        local z = cell(nargs)
        for i = 1:nargs
            z[i] = expr(:call, { :*, dargs[i], args[1:i-1]..., args[i+1:end]...})
        end
        expr(:call, { :+, z... })
    elseif op == :-
        @assert nargs <= 2
        if nargs == 1
            :( - $(dargs[1]) )
        else
            :( $(dargs[1]) - $(dargs[2]) )
        end
    elseif op == :/
        @assert nargs == 2
        :( ($(dargs[1])*$(args[2]) - $(args[1])*$(dargs[2]))/square($(args[2])) )
    elseif op == :exp
        @assert nargs == 1
        :( $(dargs[1])*exp($(args[1])) )
    elseif op == :log
        @assert nargs == 1
        :( $(dargs[1])/$(args[1]) )
    elseif op == :square
        @assert nargs == 1
        :( 2*$(dargs[1])*$(args[1]) )
    elseif op == :^
        @assert nargs == 2
        if dargs[2] == 0
            return :( $(dargs[1])*$(args[2])*$(args[1])^($(args[2])-1.0) )
        else
            error("Unimplemented")
        end
    else
        error("Unhandled operator $(op)")
    end
end


# Replace an expression a=b by a-b
homogeneize(e::Number) = e
homogeneize(e::Symbol) = e
function homogeneize(e::Expr)
    if e.head == :(=)
        expr(:call, { :-, e.args... })
    else
        e
    end
end

# Converts a dynamic expression to its static equivalent
to_static(e::Number, dynvars::Vector{Symbol}) = e
to_static(e::Symbol, dynvars::Vector{Symbol}) = e
function to_static(e::Expr, dynvars::Vector{Symbol})
    @assert e.head == :call
    if contains(dynvars, e.args[1])
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        e.args[1]
    else
        expr(:call, { e.args[1], map(x -> to_static(x, dynvars), e.args[2:end])... })
    end
end

# Substitutes a symbol by some expression
subst_symb(e::Number, env::Dict{Symbol, MathExpr}) = e
subst_symb(e::Symbol, env::Dict{Symbol, MathExpr}) = has(env, e) ? env[e] : e
function subst_symb(e::Expr, env::Dict{Symbol, MathExpr})
    @assert e.head == :call
    expr(:call, { e.args[1], map(x -> subst_symb(x, env), e.args[2:end])... })
end
 
# Substitutes leads and lags by symbols (with no lead/lag) given in dictionaries
subst_lag_lead(e::Number, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol}) = e
subst_lag_lead(e::Symbol, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol}) = e
function subst_lag_lead(e::Expr, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol})
    @assert e.head == :call
    if has(lag_dict, e.args[1]) && e.args[2] == -1
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        lag_dict[e.args[1]]
    elseif has(lead_dict, e.args[1]) && e.args[2] == 1
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        lead_dict[e.args[1]]
    else
        expr(:call, { e.args[1], map(x -> subst_lag_lead(x, lag_dict, lead_dict), e.args[2:end])... })
    end
end

