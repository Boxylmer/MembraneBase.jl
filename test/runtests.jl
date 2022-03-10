using MembraneBase
using Test

using Measurements


@testset "MembraneBase.jl" begin
    
    # ensure_matrices_are_same_size
    mat1 = rand(1, 2, 3)
    mat2 = rand(1, 2, 3)
    mat3 = rand(1, 3, 3)
    @test ensure_matrices_are_same_size(mat1, mat2) == size(mat1)
    @test_throws DimensionMismatch ensure_matrices_are_same_size(mat1, mat3)


    # add_uncertainty_to_values
    val = 4
    @test add_uncertainty_to_values([val], 0.1)[1] == 4 ± 0.4

    # chop_vector_at_first_missing_value
    myvec = [1, 2, missing, 3]
    @test [1, 2] == chop_vector_at_first_missing_value(myvec)

    # first_nonmissing_parameter
    myvec = [missing, missing, 2, missing, 3]
    @test first_nonmissing_parameter(myvec...) == 2
    @test ismissing(first_nonmissing_parameter(missing, missing))

    # safe_parse_bool
    @test safe_parse_bool.(["True", "False", "true", "false", "t, ", "gf", missing]) == [true, false, true, false, false, false, false]

    # strip_measurement_to_value
    vals = [nothing, 1, 1 ± 3, [1 ± 1, 1 ± 1.12], [1 1; 1 1 ± 1]]
    @test strip_measurement_to_value.(vals) == [nothing, 1, 1, [1, 1], [1 1; 1 1]]

    # contains_measurement_type
    @test contains_measurement_type([1, 2, 3 ± 2])
    
    # squared_uncertainty_contribution
    # fractional_uncertainty_contribution
    a = 1 ± 2
    b = 1 ± 4
    c = a + b
    @test squared_uncertainty_contribution(c, b) == 16
    @test round(fractional_uncertainty_contribution(c, a); digits=8) == 0.2

    # rss
    @test rss([2, 3], [0.8, 3]) == (2-0.8)^2
    @test rss([2±2, 3±2], [0.8±1, 3]) == 0.36

    # approximate_hessian
    # inverse_hessian
    # rss_standard_error_squared
    # rss_covariance_matrix
    # rss_minimizer_standard_errors
        # TODO 

    
    # fitting linear data
    x = [1, 2, 3]
    y = [2, 4.1, 5.8]
    y_meas = [2 ± 0.0, 4.1 ± 0.4, 5.8 ± 0.3]
    @test_logs (:warn,) fit_linear_data(x, y_meas)  # using a measurement with perfect precision should throw a warning
    @test_nowarn fit_linear_data(x, y)  

    # conversion methods
    molecular_weights = [1, 3., 4.13, π, 1e30]
    mole_fractions = [0.2, 0.4, 0.1, 0.2, 0.1]

    mass_fractions = mole_fractions_to_mass_fractions(mole_fractions, molecular_weights)
    recovered_mole_fractions = mass_fractions_to_mole_fractions(mass_fractions, molecular_weights)
    @test mole_fractions ≈ recovered_mole_fractions

    molar_volume = 0.228  # L/mol
    density = molar_volume_to_density(molar_volume, mole_fractions, molecular_weights)
    recovered_molar_volume = density_to_molar_volume(density, mole_fractions, molecular_weights)
    @test molar_volume ≈ recovered_molar_volume

    # taylor series expansions
    function myfunc(x) 
        return [
            log(x[1] + x[2] + x[3]+ x[4]+ x[5] + 0.5),
            log(x[1] + x[2] + x[3]+ x[4]+ x[5] + 0.4),
            log(x[1] + x[2] + x[3]+ x[4]+ x[5] + 0.3),
            log(x[1] + x[2] + x[3]+ x[4]+ x[5] + 0.2),
            log(x[1] + x[2] + x[3]+ x[4]+ x[5] + 0.1),
            ]
    end
    taylorfunc = get_taylor_series_function(myfunc, [0.5, 0.5, 0.5, 0.5, 0.5], 17)
    @test taylorfunc([0.1, 0.5, 0.1, 0.2, 0.56]) ≈ myfunc([0.1, 0.5, 0.1, 0.2, 0.56])
   


end
