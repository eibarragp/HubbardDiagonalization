module Utils

"""
    convert_strings_to_symbols(dict::Dict{String,Any}) -> Dict{Symbol,Any}

Convert the keys of a dictionary from `String` to `Symbol`. This allows a dictionary loaded from a TOML file to be used as a keyword argument list.
"""
function convert_strings_to_symbols(dict::Dict{String,Any})
    new_dict = Dict{Symbol,Any}()
    for (key, value) in dict
        new_dict[Symbol(key)] = value
    end
    return new_dict
end

"""
    map_dict_values(f::Function, dict::Dict{K,V}) -> Dict{K,Any}

Applies a function to the values of a dictionary, returning a new dictionary with the same keys and the transformed values.
f: The function to apply to the values.
dict: The dictionary whose values should be transformed.
"""
map_dict_values(::Type{T}, f::Function, dict::Dict{K,V}) where {K,V,T} =
    Dict{K,T}((k, T(f(v))) for (k, v) in dict)
map_dict_values(f::Function, dict::Dict{K,V}) where {K,V} = map_dict_values(Any, f, dict)

end
