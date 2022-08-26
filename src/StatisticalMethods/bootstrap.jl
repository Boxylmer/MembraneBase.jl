# realy just a wrapper module to smooth the transition to the foreign bootstrap library
function bootstrap_uncertainty_original(fitting_function, data; nsamples=500)
    results = bootstrap(
        fitting_function, 
        data, 
        BalancedSampling(nsamples)
    )
    return stderror(results)
end

evenly_sample(vec, n; replacement=true) = @views sample(vec[:], n; replace=replacement)
evenly_sample!(preallocated, vec; replacement=true) = @views sample!(vec[:], preallocated; replace=replacement)

function bootstrap_uncertainty(fitting_function, data::AbstractVector{T}; nsamples=500, replacement=true) where {T}
    # fittings function(vec) -> results
    base_results = fitting_function(data)
    # prealloc_indices = MVector{length(data), Int64}(undef)
    prealloc_inputs = MVector{length(data), T}(undef)
    variances = SVector{length(base_results), Variance}(Variance() for _ in eachindex(base_results))
    prealloc_results = MVector{length(base_results), Float64}(undef)

    for idx in 1:nsamples
        evenly_sample!(prealloc_inputs, data; replacement)
        prealloc_results .= fitting_function(prealloc_inputs)
        for idx in eachindex(prealloc_results, variances)
            fit!(variances[idx], prealloc_results[idx])
        end
        # prealloc_results[idx] = SVector(fitting_function(prealloc_inputs))
    end

    return std.(variances)
    # stdev = std.(zip(prealloc_results...))
    # return stdev
    # f(x,n; replace) = @views sample(CartesianIndices(x[:]), n; replace)

    # return stderror(results)
end


# function fitfunc(x::AbstractVector)

#     return (1 + randn(1)[1], 1 + randn(1)[1] * 5, 1 + randn(1)[1] * 9)
# end

# bootstrap_uncertainty(fitfunc, [(.1, 2.), (2., 4.)]; nsamples=5000)
# @btime bootstrap_uncertainty($fitfunc, rand(1000); nsamples=1000)