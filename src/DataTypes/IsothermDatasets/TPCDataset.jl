struct TPCDataset{TPCVT} <: AbstractArray{eltype(TPCVT), 1}
    tpc_vectors::TPCVT
    using_fugacity::Bool
    # temperatures::TVT 
end

function TPCDataset(isotherms::AbstractVector{<:IsothermData}; use_fugacity::Bool = false)
    for isotherm in isotherms
        if !(num_components(isotherm) == 1)
            throw(ErrorException("Not implemented yet for multicomponent isotherms. Components: " * string(num_components(isotherm))))
        end
    end
    # total_steps = sum(SVector{length(isotherms), Int64}[num_steps(isotherm) for isotherm in isotherms])

    # tpc_type = promote([typeof(temperature(isotherm))])
    # tpc_vecs = SVector{total_steps, Vector}()
    # for isotherm in isotherms

    # end
    tpc_vecs = []
    for isotherm in isotherms
        # determined staticarrays are about 3x faster for allocation time
        # todo try using svectors on these
        if use_fugacity
            pressures = fugacities(isotherm; component=1)
        else
            pressures = partial_pressures(isotherm; component=1)
        end
        append!(tpc_vecs, [@SVector[temperature(isotherm), pressures[idx], concentration(isotherm; component=1, step=idx)] for idx in 1:num_steps(isotherm)])
    end
    sort!(tpc_vecs)
    return TPCDataset(tpc_vecs, use_fugacity)

    # if we make the vector [[T1, P1, C1], [T2, P2, C2], ...] etc, then it's difficult to make P and C vectors as well and still do things with staticarrays
end

TPCDataset(isotherm::IsothermData; use_fugacity::Bool = false) = TPCDataset([isotherm]; use_fugacity)

function get_isotherms(tpc_dataset::TPCDataset)
    function add_tpc_vecs_to_isotherm_array!(isotherm_vec::AbstractVector{<:IsothermData}, const_t_tpc_vecs::AbstractVector, using_fugacity::Bool)
        t_vec, p_vec, c_vec = zip(const_t_tpc_vecs...)
        if using_fugacity
            iso = IsothermData(; concentrations_cc = c_vec, fugacities_mpa = p_vec, temperature_k = t_vec[1])
        else
            iso = IsothermData(; concentrations_cc = c_vec, partial_pressures_mpa = p_vec, temperature_k = t_vec[1])
        end
        push!(isotherm_vec, iso)
    end
    
    isotherms = IsothermData[]
    current_temperature_index = 1
    current_temperature = tpc_dataset.tpc_vectors[current_temperature_index][1]
    num_tpc_vecs = length(tpc_dataset.tpc_vectors)

    for (tpc_vec_idx, tpc_vec) in enumerate(tpc_dataset)
        # Check if we're at the last index
            # if we are,
                # check if [current index] is the current temperature
                    # if it is, slice from the last found index to the end
                    # if it isn't make a singleton isotherm
            # if we are not, then 
                # check if [current index + 1] is the same temperature
                    # if it is, continue
                    # if not, slice at [current_temperature_index:tpc_vec_idx] and make the isotherm
                        # also set the current temperature to the index + 1's temperature

        if tpc_vec_idx != num_tpc_vecs  # if we're not at the end of the iteration
            if tpc_dataset.tpc_vectors[tpc_vec_idx + 1][1] == current_temperature  # if the next item in the iteration is the same temperature
                continue
            else # if the next item in the iteration is a different temperature
                vecs_slice = tpc_dataset.tpc_vectors[current_temperature_index:tpc_vec_idx]
                current_temperature_index = tpc_vec_idx + 1
                current_temperature = tpc_dataset.tpc_vectors[current_temperature_index][1]
            end
        else  # if we're at the end of the iteration
            if current_temperature == tpc_vec[1] # if the last item is the same temerature as before
                vecs_slice = tpc_dataset.tpc_vectors[current_temperature_index:tpc_vec_idx]
            else  # if the last item is a unique temperature
                vecs_slice = tpc_dataset.tpc_vectors[tpc_vec_idx:tpc_vec_idx]
            end
        end

        add_tpc_vecs_to_isotherm_array!(isotherms, vecs_slice, tpc_dataset.using_fugacity)

    end
    return isotherms
end

# AbstractArray implementation
# tpc_vectors looks like [[T1, P1, C1], [T2, P2, C2], ...] where T, P, and C, are some <:Number
function Base.size(obj::TPCDataset)
    return size(obj.tpc_vectors)
end

function Base.getindex(obj::TPCDataset, i::Int) 
    return obj.tpc_vectors[i]
end

Base.length(obj::TPCDataset) = length(obj.tpc_vectors)

Base.copy(obj::TPCDataset) = TPCDataset(copy(obj.tpc_vectors), copy(obj.using_fugacity))
Base.deleteat!(obj::TPCDataset, index) = deleteat!(obj.tpc_vectors, index)
