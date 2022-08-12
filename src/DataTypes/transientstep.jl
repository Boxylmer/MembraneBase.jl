struct TransientStepData{N}
    time::MVector{N, Float64}  # seconds
    dimensionlesssorption::MVector{N, Float64}  # unitless
end
function TransientStepData(t::AbstractVector, dimensionlesssorption::AbstractVector)
    n = length(t)
    @assert n == length(dimensionlesssorption)
    tt = eltype(t)
    st = eltype(dimensionlesssorption)
    if st <: Measurement && tt <: Measurement
        return TransientStepData(
            MVector{n, Float64}(strip_measurement_to_value(t)), 
            MVector{n, Float64}(strip_measurement_to_value(dimensionlesssorption))
        )
    elseif st <: Measurement
        return TransientStepData(MVector{n, Float64}(t), MVector{n, Float64}(strip_measurement_to_value(dimensionlesssorption)))
    elseif tt <: Measurement
        return TransientStepData(MVector{n, Float64}(strip_measurement_to_value(t)), MVector{n, Float64}(dimensionlesssorption))
    else
        return TransientStepData(MVector{n, Float64}(t), MVector{n, Float64}(dimensionlesssorption))
    end
end



"""
    TransientStepData(dataset)
Container for data concerning a single transient sorption step. I.e., a profile of mass concentration with time (starting from a previous equilibrium and *generally* reaching a new equilibrium) of some penetrant into an even film of polymer material, or more succinctly, a *polymer slab*. 

# Arguments
- `dataset`: Should be some iterable with 
    - the first item being an iterable of time data (in **seconds**)
    - the second item being an iterable of dimensionless sorption data (no units)  
"""
function TransientStepData(dataset::Union{Tuple, AbstractVector})
    if eltype(dataset) <: Tuple{Float64, Float64}
        reinterpretation = reinterpret(reshape, Float64, vec)
        time_data = @view reinterpretation[1, :]
        sorption_data = @view reinterpretation[:, 1]
    else
        @warn "Non-Float64 types handed to TransientStepData, conversion will be very slow."
        time_data = [datum[1] for datum in dataset]
        sorption_data = [datum[2] for datum in dataset]
    end
    return TransientStepData(time_data, sorption_data)
end

"""
    resample(transient::TransientStepData, num_datapoints, time_function)
Resample the data with some `time_function` spacing and linear interpolation between adjacent data points.
Use this method when you have a ridiculous amount of data and need to slim it down while losing as little valuable information as possible.

# Arguments
- `transient::TransientStepData`: TransientStepData used to resample
- `num_datapoints::Integer`: number of resampled data points to be in the resulting TransientStepData (should always be less than what you currently have)
- `time_function::Symbol`: how to weigh the time data...
    * `:Linear` will have even spacing from t=0 to t=num_datapoints
    * `:Root` will sample along the root of time (more information about the initial behavior will be kept)
        * It is usually best to use :Root for most scenarios and :Linear for quick scenarios
    * `:Log` will sample along the log_10 of time (much more info about initial behavior will be kept)
        * Reserve Log for extremely long timescales (think: on the order of *months* or *years*, not *days* or *weeks*)

This function assumes that time data is sorted in ascending order
"""
function resample(transient::TransientStepData, num_datapoints, time_function; skip_validity_checks=false)
    # todo: remove any data point at or below zero if using :Log, right now theres a quick fix that just removes 0
    min_time, max_time = extrema(transient.time)
    
    if !skip_validity_checks
        if minimum(transient.time[i+1] - transient.time[i] for i in eachindex(transient.time)[1:end-1])[1] < 0
            idx = findmin(collect(transient.time[i+1] - transient.time[i] for i in eachindex(transient.time)[1:end-1]))[2]
            throw(DomainError(
                "Warning: Time must always increase. At least one point in the transient step was backwards in time (largest change at index " * string(idx) * ")"))
        end
    end

    if time_function == :Linear
        # it's already linear!
        step_size = (max_time - min_time) / (num_datapoints - 1)  # off by one errors are the worst
        resampled_time = collect(min_time:step_size:max_time)
    elseif time_function == :Root
        max_time = sqrt(max_time)
        min_time = sqrt(min_time)
        step_size = (max_time - min_time) / (num_datapoints - 1)  # off by one errors are the worst
        resampled_time = collect(min_time:step_size:max_time).^2
    elseif time_function == :Log
        # remove zero if present
        if min_time == 0
            min_time = minimum(transient.time[2:end])
        end
        max_time = log10(max_time)
        min_time = log10(min_time)
        step_size = (max_time - min_time) / (num_datapoints - 1)
        resampled_time = 10 .^ (collect(min_time:step_size:max_time))
    elseif isnothing(time_function)
        return transient
    else
        throw(ArgumentError("A valid time function was not given (try \":Root\"?)"))
    end
    # some weirdness with floating point errors can occur here that makes searching by index a pain
    # here, we make absolutely sure the values are *exactly* what they should be
    resampled_time[1] = transient.time[1]
    resampled_time[end] = transient.time[end]     
    
    # resampled_sorption = similar(resampled_time, eltype(transient.dimensionlesssorption))
    resampled_sorption = similar(resampled_time)
    
    for (sorption_idx, t) in enumerate(resampled_time)
        idx1 = searchsortedfirst(transient.time, t)
        if transient.time[idx1] == t
            resampled_sorption[sorption_idx] = transient.dimensionlesssorption[idx1]
        else
            idx0 = idx1 - 1
            interpolated_sorption = interpolate_linearly(transient.time[idx0], transient.dimensionlesssorption[idx0], transient.time[idx1], transient.dimensionlesssorption[idx1], t)
            resampled_sorption[sorption_idx] = interpolated_sorption
        end
    end
    return TransientStepData(resampled_time, resampled_sorption)

end

function dataset(transient::TransientStepData)
    
    return collect(zip(transient.time, transient.dimensionlesssorption))
end

function strip_measurement_to_value(meas::TransientStepData)
    return TransientStepData(strip_measurement_to_value(meas.time), strip_measurement_to_value(meas.dimensionlesssorption))
end