# realy just a wrapper module to smooth the transition to the foreign bootstrap library
function bootstrap_uncertainty(fitting_function, data; nsamples=500)
    results = bootstrap(
        fitting_function, 
        data, 
        BalancedSampling(nsamples)
    )
    return stderror(results)
end

# evenly_sample(vec, n; replacement=true) = @views sample(vec[:], n; replace=replacement)
# evenly_sample!(preallocated, vec; replacement=true) = @views sample!(vec[:], preallocated; replace=replacement)

# function bootstrap_uncertainty(fitting_function, data::AbstractVector{T}; nsamples=500, replacement=true) where {T}
#     # fittings function(vec) -> results
#     base_results = fitting_function(data)

#     # prealloc_indices = MVector{length(data), Int64}(undef)
#     prealloc_inputs = MVector{length(data), T}(undef)
#     prealloc_results = MVector{nsamples, SVector{length(base_results), Float64}}(undef)

#     for idx in 1:nsamples
#         evenly_sample!(prealloc_inputs, data; replacement)
#         prealloc_results[idx] = SVector(fitting_function(prealloc_inputs))
#     end
    
#     stdev = std.(zip(prealloc_results...))
#     println(prealloc_results)
#     # f(x,n; replace) = @views sample(CartesianIndices(x[:]), n; replace)

#     # return stderror(results)
# end


# function fitfunc(x::AbstractVector)

#     return [1 + randn(1)[1], 1 + randn(1)[1] * 5]
# end