# Replace an expression a=b by a-b
homogeneize(e::Number) = e
homogeneize(e::Symbol) = e
function homogeneize(e::Expr)
    if e.head == :(=)
        Expr(:call, :-, e.args...)
    else
        e
    end
end

# Converts a dynamic expression to its static equivalent
to_static(e::Number, dynvars::Vector{Symbol}) = e
to_static(e::Symbol, dynvars::Vector{Symbol}) = e
function to_static(e::Expr, dynvars::Vector{Symbol})
    @assert e.head == :call
    if e.args[1] in dynvars
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        e.args[1]
    else
        Expr(:call, e.args[1], map(x -> to_static(x, dynvars), e.args[2:end])...)
    end
end

# Substitutes a symbol by some expression
subst_symb(e::Number, env::Dict{Symbol, MathExpr}) = e
subst_symb(e::Symbol, env::Dict{Symbol, MathExpr}) = haskey(env, e) ? env[e] : e
function subst_symb(e::Expr, env::Dict{Symbol, MathExpr})
    @assert e.head == :call
    Expr(:call, e.args[1], map(x -> subst_symb(x, env), e.args[2:end])...)
end
 
# Substitutes leads and lags by symbols (with no lead/lag) given in dictionaries
subst_lag_lead(e::Number, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol}) = e
subst_lag_lead(e::Symbol, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol}) = e
function subst_lag_lead(e::Expr, lag_dict::Dict{Symbol,Symbol}, lead_dict::Dict{Symbol,Symbol})
    @assert e.head == :call
    if haskey(lag_dict, e.args[1]) && e.args[2] == -1
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        lag_dict[e.args[1]]
    elseif haskey(lead_dict, e.args[1]) && e.args[2] == 1
        @assert length(e.args) == 2
        @assert isa(e.args[2], Integer)
        lead_dict[e.args[1]]
    else
        Expr(:call, e.args[1], map(x -> subst_lag_lead(x, lag_dict, lead_dict), e.args[2:end])...)
    end
end

