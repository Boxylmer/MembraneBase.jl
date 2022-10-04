module MembraneBase
    using Measurements
    using FiniteDiff
    using TaylorSeries
    using GLM
    using DataFrames
    using Bootstrap
    using StaticArrays
    using LinearAlgebra
    using StatsBase: sample!, sample
    # using Statistics
    using OnlineStats
    
    include("Constants.jl")
    include("HelperFunctions.jl")

    # main utility methods
    export ensure_matrices_are_same_size
    export add_uncertainty_to_values
    export chop_vector_at_first_missing_value
    export first_nonmissing_parameter
    export safe_parse_bool
    export strip_measurement_to_value
    export contains_measurement_type
    export squared_uncertainty_contribution
    export fractional_uncertainty_contribution
    export rss
    export approximate_hessian
    export inverse_hessian
    export rss_standard_error_squared
    export rss_covariance_matrix
    export rss_minimizer_standard_errors
    export estimate_slope_by_adjacent_points
    export interpolate_linearly
    export fit_linear_data
    export checkreal
    export filterreal
    export solve_cubic_eq
    export get_taylor_series_function
    export expand_taylor_higher_terms
    export antoines_vapor_pressure  # todo deprecate this
    export mass_fractions_to_mole_fractions
    export mole_fractions_to_mass_fractions
    export molar_volume_to_density
    export density_to_molar_volume
    export polymer_phase_mass_fractions_to_gpen_per_gpol
    export polymer_phase_mass_fractions_to_ccpen_per_ccpol
    export polymer_specific_volume

    # isotherm methods
    include(joinpath("DataTypes", "isotherm.jl"))
    export IsothermData
    export isotherm_dataset
    export pressure
    export partial_pressures
    export concentration
    export mole_fractions
    export activities
    export fugacities
    export temperature
    export polymer_density
    export molecular_weights
    export num_components
    export num_steps
    export mass_sorbed
    export penetrant_mass_fractions
    export increasing_concentration
    export remove_increasing_concentration_steps

    # transient step data
    include(joinpath("DataTypes", "transientstep.jl"))
    export TransientStepData
    export resample
    export dataset
    
    # statistical methods
    include(joinpath("StatisticalMethods", "bootstrap.jl"))
    include(joinpath("StatisticalMethods", "jackknife.jl"))
    export bootstrap_uncertainty
    export jackknife_uncertainty

end
