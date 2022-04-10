# realy just a wrapper module to smooth the transition to the foreign bootstrap library
function bootstrap_uncertainty(fitting_function, data; nsamples=500)
    results = bootstrap(
        fitting_function, 
        data, 
        BalancedSampling(nsamples)
    )
    return stderror(results)
end