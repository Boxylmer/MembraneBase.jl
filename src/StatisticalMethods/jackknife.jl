
"""
    jackknife_uncertainty(fitting_function, data)

    Requirements:
        `data` can be any iterable, including custom structs.
        `fitting_function` must accept the the same iterable (data) which will be resampled, 
            it must return a vector of the fitted parameters.
    This method will return a vector of uncertainties in the same order as the fitted parameters.
"""
function jackknife_uncertainty(fitting_function::Function, data::AbstractVector)
    resampled_fittings = _jackknife_resampling_and_fitting(fitting_function, data)
    true_fitting = fitting_function(data)
    variances = [_jackknife_variance(resampled_fitting, true_fitting[idx]) for (idx, resampled_fitting) in enumerate(resampled_fittings)]
    uncertainties = variances .^(0.5)
    return uncertainties
end

function jackknife_bias(fitting_function::Function, data::AbstractVector)
    num_datapoints = length(data)
    resampled_fittings = _jackknife_resampling_and_fitting(fitting_function, data)
    true_fitting = fitting_function(data)
    average_resampled_fitting = [sum(fittings) for fittings in resampled_fittings] ./ num_datapoints
    
    biases = (num_datapoints - 1) .* (average_resampled_fitting - true_fitting)
    return biases
end

function _jackknife_resampling_and_fitting(fitting_function::Function, data::AbstractVector)
    num_datapoints = length(data)
    resampled_fittings = Vector{Vector{Any}}(undef, num_datapoints)
    for counter_id in 1:num_datapoints
        resampled_data = _resample(data, counter_id)
        resampled_output = fitting_function(resampled_data)
        resampled_fittings[counter_id] = resampled_output
    end
    return collect(zip(resampled_fittings...))
end

function _resample(data::AbstractVector, idx::Integer)
    new_data = copy(data)
    deleteat!(new_data, idx)
    return new_data
end

function _jackknife_variance(resampled_values, real_value)
    num_pts = length(resampled_values)
    return sum((resampled_values .- real_value).^2) * (1/(num_pts - 1))
end