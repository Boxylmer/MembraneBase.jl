
# utility methods

"Throw an error if the given matrixes aren't the same size, otherwise return their size."
function ensure_matrices_are_same_size(mats...)
    discovered_size = nothing
    for mat in mats  # for every matrix to be compared,
        if !isnothing(mat)  # if there is data in the matrix, 
            if !isnothing(discovered_size)  # and if we already know the size we're comparing,
                if !(discovered_size == size(mat))  # do the comparison and throw an error if theres a mismatch,
                    throw(DimensionMismatch("Dimensions of steps and components did not match in one or more inputs."))
                end
            else  # otherwise, we should use this size as the comparison in the future
                discovered_size = size(mat)
            end 
        end
    end
    return discovered_size  # hand the verified size back for convenience
end

"add an uncertainty of m*x to the values x where m is the uncertainty_multiplier."
function add_uncertainty_to_values(values::AbstractVector, uncertainty_multiplier::Number)
    return [values[i] ± (values[i] * uncertainty_multiplier) for i in eachindex(values)]
end

"Run through vector and cut it at the first missing value."
function chop_vector_at_first_missing_value(vector::AbstractVector)
    for i in eachindex(vector)
        if ismissing(vector[i])
            return vector[1:i-1]
        end
    end
end

"Return the first value that isn't missing, if all are missing, return missing."
function first_nonmissing_parameter(args...)
    for arg in args
        if !ismissing(arg)
            return arg
        end
    end
    return missing
end

"Safe parse boolean."
function safe_parse_bool(string_to_parse::AbstractString; default=false) 
    formatted_string = strip(lowercase(string_to_parse), [' ',])
    if formatted_string == "true" || formatted_string == "t"
        return true
    else
        return default
    end
end
safe_parse_bool(::Missing; kwargs...) = safe_parse_bool(" "; kwargs...)
safe_parse_bool(input::Bool; kwargs...) = input

# handling measurement types
"Return the base value of the measurement, regardless of what kind of measurement it is. If it isn't a measurement, return the object passed. If it's an array or struct containing measurements, return an equivalent array or struct of values."
function strip_measurement_to_value(::Nothing)
    return nothing
end
function strip_measurement_to_value(::Missing)
    return missing
end
function strip_measurement_to_value(obj::Measurement)
    measurement_value = obj.val
    return measurement_value
end
function strip_measurement_to_value(obj::AbstractFloat)
    return obj
end
function strip_measurement_to_value(obj::Integer)
    return obj
end
function strip_measurement_to_value(vector::AbstractVector{<:Measurement})
    return [vector[i].val for i in eachindex(vector)]
end
function strip_measurement_to_value(measurement_matrix::AbstractMatrix{<:Measurement})
    value_matrix = map(x -> x.val, measurement_matrix)
    return value_matrix
end
function strip_measurement_to_value(iterable::AbstractArray{<:Number})
    return iterable
end
function strip_measurement_to_value(weird_iterable::AbstractArray{<:Any})
    return strip_measurement_to_value.(weird_iterable)
end

"Check if any of the arguments contain measurements."
function contains_measurement_type(args...)
    for arg in args
        if contains_measurement_type(arg)
            return true
        end
    end
    return false
end
function contains_measurement_type(arg::Number)
    if typeof(arg) <: Measurement return true end
    return false
end
function contains_measurement_type(arg::AbstractArray)
    if eltype(arg) <: Measurement
        return true
    end
    return false
end

"Calculate the (squared) contribution a measurement `var` has on a `result`ing measurement."
function squared_uncertainty_contribution(result::Measurement, var::Measurement)
    var_key = (var.val, var.err, var.tag)

    components = Measurements.uncertainty_components(result)
    if haskey(components, var_key)
        drdv_times_sigvar = components[var_key]
        drdv_times_sigvar_squared = drdv_times_sigvar^2
        return drdv_times_sigvar_squared
    else
        return 0.0
    end
end
function squared_uncertainty_contribution(::Measurement, ::Number) 
    return 0.0 # if the variable was not a measurement, the contribution to error is always zero
end

"Calculate the fractional contribution a measurement `var` has on a `result`ing measurement."
function fractional_uncertainty_contribution(result::Measurement, var::Number)
    return squared_uncertainty_contribution(result, var)/result.err^2 
end

# approximations and general mathy things

"Calculate the sum of squared residuals of an `experimental` value to some `predicted` value`."
function rss(experimental::AbstractVector{<:Measurement}, predicted::AbstractVector{<:Number})
    return sum(([item.val for item in experimental] - strip_measurement_to_value(predicted)).^2 ./ [item.err for item in experimental].^2)
end
function rss(experimental::AbstractVector{<:Number}, predicted::AbstractVector{<:Number})
    return sum((strip_measurement_to_value(experimental) - strip_measurement_to_value(predicted)).^2)
end

# Find the parameter error and other things involving RSS 
"Estimate the hessian matrix of a function with respect to a vector of values. The value vector should match the input of the `func`tion."
function approximate_hessian(func::Function, values::AbstractVector{<:Number})
    return FiniteDiff.finite_difference_hessian(func, values)
end

# todo as of 8/25/2021, symmetric staticarrays do not work with inv(). Thus, we need to return it to a normal matrix for now. If this issue is fixed, remove this and test.
# https://github.com/JuliaArrays/StaticArrays.jl/issues/955 
"Calculate the inverse hessian matrix of an rss function `rss_func`, with respect to some minimizer `minimizer_values`."
inverse_hessian(rss_func::Function, minimizer_values::AbstractVector{<:Number}) = inv(Symmetric(Matrix(approximate_hessian(rss_func, minimizer_values))))

"Estimate the standard error squared" # depricated? todo remove this
rss_standard_error_squared(rss_value, num_params::Integer, num_datapts::Integer) = (rss_value / (num_datapts - num_params))

# function rss_covariance_matrix(rss_func::Function, minimizer_values::AbstractVector{<:Number}, num_observations::Integer)
#     inv_hess = inverse_hessian(rss_func, minimizer_values)
#     rss_squared = rss_standard_error_squared(rss_func(minimizer_values), length(minimizer_values), num_observations)
#     cov_mat = inv_hess * rss_squared
#     return cov_mat
# end 

"Estimate the RSS covariance matrix using the inverse hessian method."
function rss_covariance_matrix(rss_func::Function, minimizer_values::AbstractVector{<:Number}, num_observations::Integer)
    inv_hess = inverse_hessian(rss_func, minimizer_values)
    # rss_squared = rss_standard_error_squared(rss_func(minimizer_values), length(minimizer_values), num_observations)
    dof = num_observations - length(minimizer_values)
    scaling_factor = rss_func(minimizer_values) / dof
    cov_mat = inv_hess * scaling_factor
    return cov_mat
end 

"Get standard errors from an RSS function using the already optimized (minimizer) values."
function rss_minimizer_standard_errors(rss_func::Function, minimizer_values::AbstractVector{<:Number}, num_observations::Integer)
    mat = rss_covariance_matrix(rss_func, minimizer_values, num_observations)
    covariance_diagonals = diag(mat)
    std_err_vec = sqrt.(abs.(covariance_diagonals))
    return std_err_vec
end

# linear fittings
"Using two adjacent points in a vector, estimate the slope of a series."
function estimate_slope_by_adjacent_points(x_vector::AbstractVector, y_vector::AbstractVector, point_index::Integer)
    lenx = length(x_vector)
    if point_index == 1
        slope = (y_vector[2] - y_vector[1])/(x_vector[2] - x_vector[1])
    elseif point_index == lenx
        slope = (y_vector[end] - y_vector[end-1])/(x_vector[end] - x_vector[end-1])
    else
        slope = (y_vector[point_index + 1] - y_vector[point_index - 1])/(x_vector[point_index + 1] - x_vector[point_index - 1])
    end
    return slope
end

function estimate_slope_by_adjacent_points(x_vector::AbstractVector, y_vector::AbstractVector)
    return [estimate_slope_by_adjacent_points(x_vector, y_vector, idx) for idx in eachindex(y_vector)]
end

"Apply linear interpolation to two points."
function interpolate_linearly(x0, y0, x1, y1, x)
    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    return y
end

"""
    fit_linear_data(x_data, y_data)

Fit two vectors of data to a linear model and return their slope and intercept. 
Both the slope and intercept will be `measurement` types, with the uncertainty being the standard error of the fittings. 

!!! note
    If the input data are `measurement`s, the `y_data` uncertainties will be used to weight the importance of each data point.
    [Alexander Aitken](https://en.wikipedia.org/wiki/Alexander_Aitken) showed that weighting according to the inverse of variance is the [best linear unbiased estimator (BLUE)](https://en.wikipedia.org/wiki/Gauss%E2%80%93Markov_theorem).
    Read more on weighted least squares in linear regression [here](https://en.wikipedia.org/wiki/Weighted_least_squares).

"""
function fit_linear_data(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}; y_weights=missing)
    data = DataFrame(X=strip_measurement_to_value(x), Y=y)
    if ismissing(y_weights)
        y_weights = fill!(similar(y), 1)
    end
    ols = lm(@formula(Y ~ X), data;  wts = y_weights)
    intercept, slope = coef(ols) .± stderror(ols)
    return slope, intercept
end

function fit_linear_data(x::AbstractVector{<:Number}, y::AbstractVector{<:Measurement})
    weights = [1 / (meas.err^2) for meas in y]
    for wt in weights
        if isinf(wt)
            @warn("A measurement uncertainty was zero. (weight_i = 1/var_i^2 = Inf when var_i = 0)\nDefaulting to unweighted linear regression.")
            weights = missing
        end
    end
    y_vals = [meas.val for meas in y]
    if eltype(x) <: Measurement x_vals = [meas.val for meas in x]
    else x_vals = x
    end
    
    return fit_linear_data(x_vals, y_vals; y_weights = weights)
end


"""
    checkreal(num::Complex; threshold = 1e-14)
Make the value real if the imaginary component is below the threshold; useful for precision errors
"""
checkreal(num::Complex; threshold = 1e-14) = abs(imag(num)) < threshold ? real(num) : num

"""
    filterreal(values::AbstractVector)
Get all values in the vector that have no imaginary component
"""
function filterreal(values::Union{AbstractVector, Tuple})
    result = []
    for val in values
        if isreal(val)
            push!(result, val)
        end
    end
    return result
end
"""
    solve_cubic_eq(poly::AbstractVector{T})
Find the roots to a cubic equation (``dx^3 + cx^2 + bx + a = 0``).
    - where poly is a vector formatted as `[a, b, c, d]`.
"""
function solve_cubic_eq(poly::AbstractVector{T}) where {T<:Real}
    # copied from PolynomialRoots.jl, adapted to be AD friendly
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    _1 = one(T)
    third = _1/3
    a1  =  complex(one(T) / poly[4])
    E1  = -poly[3]*a1
    E2  =  poly[2]*a1
    E3  = -poly[1]*a1
    s0  =  E1
    E12 =  E1*E1
    A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B   =  E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = (A*A - 4*B*B*B)^0.5
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return (
        checkreal(third*(s0 + s1 + s2)), 
        checkreal(third*(s0 + s1*zeta2 + s2*zeta1)), 
        checkreal(third*(s0 + s1*zeta1 + s2*zeta2))
        )
end

# this method is currently unused. Uncomment this if it becomes used
"""Calculate the multidimensional taylor series of a function `func` starting at the vector `x0` using `order` terms."""
function get_taylor_series_function(func::Function, x0::AbstractVector, order::Integer)
    set_variables("x", numvars=length(x0), order=order)
    expansions = taylor_expand(func, x0; order=order)
    result = function resfunc(x)
        return [expansion(x .- x0) for expansion in expansions]
    end
    return result
end

# this method is currently unused. Uncomment this if it becomes used
"""Calculate the taylor series expansion starting at `x0` of the singly valued `func` of `order` terms."""
function get_taylor_series_function(func::Function, x0::Number, order::Integer)
    result = function resfunc(x)
        return taylor_expand(func, x0; order)(x)
    end
    return result
end

"""Calculate the taylor series expansion starting at `x0` of the singly valued `func` of `order` terms without the first expansion term."""
function expand_taylor_higher_terms(func::Function, x0::Number, order::Integer)
    result = function resfunc(x::Number)
        return taylor_expand(func, x0; order)(x) - taylor_expand(func, x0; order=0)(x)
    end
    return result
end

"""Calculate the taylor series expansion starting at `x0` of the singly valued `func` of `order` terms without the first expansion term."""
function expand_taylor_higher_terms(func::Function, x0::AbstractVector, order::Integer)
    set_variables("x", numvars=length(x0), order=order)
    expansion = taylor_expand(func, x0; order=order)
    expansion_first_term = taylor_expand(func, x0; order=0)
    resulting_function = function resfunc(x::AbstractVector)
        return expansion(x .- x0) - expansion_first_term(x .- x0)
    end
    return resulting_function
end

# chemistry things (may be moved later)
function antoines_vapor_pressure(a, b, c, temp_k)
    return 10^(a - b / (c + temp_k))
end

"Convert mass fraction to mole fraction."
function mass_fractions_to_mole_fractions(mass_fractions, molecular_weights)
    moles = mass_fractions ./ molecular_weights
    mole_fractions = moles ./ sum(moles)
    return mole_fractions
end

"Convert mole fraction to mass fraction."
function mole_fractions_to_mass_fractions(mole_fractions, molecular_weights)
    masses = mole_fractions .* molecular_weights
    mass_fractions = masses ./ sum(masses)
    return mass_fractions
end

"Convert molar volume to density in g/cm^3."
function molar_volume_to_density(molar_volume_l_mol, mole_fractions, molecular_weights_g_mol)
    molar_density = 1 / molar_volume_l_mol  # L/mol -> mol/L
    mass_concentrations = molar_density .* mole_fractions .* molecular_weights_g_mol
    return sum(mass_concentrations) / 1000 # g / cm3
end

"Convert density to molar volume in L/mol."
function density_to_molar_volume(density_g_cm3, mole_fractions, molecular_weights_g_mol)
    mass_fracs = mole_fractions_to_mass_fractions(mole_fractions, molecular_weights_g_mol)
    molar_densities = density_g_cm3 * 1000 .* mass_fracs ./ molecular_weights_g_mol
    molar_volume = 1 / sum(molar_densities)
    return molar_volume  # L / mol
end

"Convert polymer phase mass fractions to `g(penetrant)/g(polymer)`. Assumes that the polymer is the first item in the vector."
function polymer_phase_mass_fractions_to_gpen_per_gpol(mass_fractions)
    return mass_fractions[2:end] ./ mass_fractions[1]
end

"Convert polymer phase mass fraction to CC(STP)/CC(polymer)."
function polymer_phase_mass_fractions_to_ccpen_per_ccpol(mass_fractions, polymer_density_g_cm3, penetrant_molecular_weights)
    # mass fractions includes the polymer
    g_per_g = polymer_phase_mass_fractions_to_gpen_per_gpol(mass_fractions)
    return g_per_g * CC_PER_MOL_STP * polymer_density_g_cm3 ./ penetrant_molecular_weights
end