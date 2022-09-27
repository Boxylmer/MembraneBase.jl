using MembraneBase
using Test
using Measurements
using BenchmarkTools
using Revise

@testset "MembraneBase.jl" begin
    

    @testset "Helper Methods" begin
        
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
        @test rss([2±2, 3±2], [0.8±1, 3]) == 1.44

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

        # the inverse hessian methods
        objective_function(xy) = ((xy[1] + xy[2])^2 + (xy[1]-1)^2 + (xy[2]-0.8)^2)^2
        minimizer = [0.494234234, 0.25034623146]
        @test errs = rss_minimizer_standard_errors(objective_function, minimizer, 10) ≈ [0.15049509481360718, 0.1517352228363271]

    end

    @testset "Statistical Methods" begin
        ndata = 100
        data_x = collect(1:ndata)
        data_y = data_x # (data_x .+ randn(100)) .^ 2.5
        data = collect(zip(data_x, data_y))
        # function predictor(b_vec)
        #     y = data_x .^ b_vec[1]
        #     err = sum((y .- data_y) .^ 2)
        #     return err / 1e9
        # end

        # this doesn't actually fit anything, it just returns 5 plus a random number
        function fitter(data_pairs)
            return [5 + randn(1)[1]]
        end

        j_σ = jackknife_uncertainty(fitter, data)
        # @test round(j_σ[1]; digits=0) == 1

        b_σ = bootstrap_uncertainty(fitter, data)
        @test round(b_σ[1]; digits=0) == 1

        # bootstrap_allocs = @allocated bootstrap_uncertainty(fitter, data; nsamples=10000)
        # @show bootstrap_allocs
        # println("Bootstrap (Will)")
        # @show @btime bootstrap_uncertainty($fitter, $data; nsamples=10000)
        # println("Bootstrap (original, Bootstrap.jl)")
        # @show @btime MembraneBase.bootstrap_uncertainty_original($fitter, $data; nsamples=10000)
    end

    @testset "Isotherm Data Structures" begin
        # Basic isotherm creastion
        pure_iso = IsothermData(partial_pressures_mpa = [1], concentrations_cc = [2])
        iso = IsothermData(partial_pressures_mpa = [[1]], concentrations_cc = [[2]])

        iso_2 = IsothermData(partial_pressures_mpa = [1, 2, 3], concentrations_cc = [2, 5, 10])
        iso_3 = IsothermData(partial_pressures_mpa = [[1, 2, 3], [1, 2, 3]], concentrations_cc = [[2, 5, 10], [1, 4, 8]])

        # invalid isotherms
        @test_throws DimensionMismatch IsothermData(partial_pressures_mpa = [1], concentrations_cc = [2, 7])
        @test_throws DimensionMismatch IsothermData(partial_pressures_mpa = [1, 2], activities = [7])
        @test_throws DimensionMismatch IsothermData(partial_pressures_mpa = [[1, 2, 3], [1, 2, 3, 4]], concentrations_cc = [[2, 5, 10, 11], [1, 4, 8]])
        @test_throws MethodError IsothermData(1, 2)
        @test_throws MethodError IsothermData(1, 2, 3)

        @test num_components(iso)   == 1
        @test num_components(iso_2) == 1
        @test num_components(iso_3) == 2

        @test pressure(iso; step = 1)   == 1
        @test pressure(iso_2; step = 3) == 3
        @test pressure(iso_3; step = 2) == 4

        # indexing behavior
        iso = IsothermData(partial_pressures_mpa = [1], concentrations_cc = [2])
        @test partial_pressures(iso; component=1)[1] == 1
        pres = [[1, 2, 3], [1.5, 2.5, 3.5]]
        conc = [[2, 5, 10], [1, 4, 8]]
        iso_2 = IsothermData(partial_pressures_mpa = pres, concentrations_cc=conc)
        @test isotherm_dataset(partial_pressures(iso_2), concentration(iso_2), 2)[2] == [2.5, 4]
        @test isotherm_dataset(partial_pressures(iso_2), concentration(iso_2))[2] == [[2.0, 2.5], [5.0, 4.0]]

        # test getter functions
        iso = IsothermData(partial_pressures_mpa = [[1, 2, 3], [1, 2, 3]], concentrations_cc = [[10, 20, 30], [40, 50, 60]], pen_mws_g_mol = [4, 6], rho_pol_g_cm3 = 1.12)
        iso_2 = IsothermData(partial_pressures_mpa = [1, 2, 3], concentrations_cc = [10, 20, 30], pen_mws_g_mol = 6, rho_pol_g_cm3 = 1.12)
        concs_g_g = concentration(iso; pol_units=:g, gas_units=:g)
        @test concs_g_g == [
            0.001593391885173807 0.009560351311042842; 
            0.003186783770347614 0.011950439138803554; 
            0.004780175655521421 0.014340526966564264]

        
        mass_fracs = penetrant_mass_fractions(iso; step=1)
        totals = [1 - sum(mass_fracs), mass_fracs...]
        g_g_converted_back =  polymer_phase_mass_fractions_to_gpen_per_gpol(totals)
        concs_g_g = concentration(iso; step=1, pol_units=:g, gas_units=:g)
        @test concs_g_g ≈ g_g_converted_back

        @test pressure(iso) == [2, 4, 6]
        @test pressure(iso; step=2) == 4
        @test pressure(iso_2) == [1, 2, 3]
        @test pressure(iso_2; step=2) == 2

        iso_3 = IsothermData(partial_pressures_mpa = [[1, 2, 3], [1, 4, 9]])
        @test mole_fractions(iso_3; step=2) == [2/6, 4/6]
        @test mole_fractions(iso_3) ==  [0.5 0.5; 2/6 4/6; 0.25 0.75]
        @test mole_fractions(iso) == [0.5 0.5; 0.5 0.5; 0.5 0.5]
        @test mole_fractions(iso; step=2) == [0.5, 0.5]
        @test mole_fractions(iso_2) == [1; 1; 1][:, :]
        

        # indexing

        # the potential issue here is that calling a getter with a singleton of either step or component takes the internal data structure, 
        #   a matrix of (steps * components) and returns a vector regardless of the original structure.
        # when reconstructing the isotherm, a vector is passed in, which is always interpreted as steps.
        pres = [[1, 2, 3, 4, 5, 6, 7], [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1], [1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2]]
        conc = [[10, 20, 30, 40, 50, 60, 70], [10.1, 20.1, 30.1, 40.1, 50.1, 60.1, 70.1], [10.2, 20.2, 30.2, 40.2, 50.2, 60.2, 70.2]]
        iso_4 = IsothermData(partial_pressures_mpa = pres, concentrations_cc = conc, pen_mws_g_mol = [100, 101, 102], rho_pol_g_cm3 = 1.12)
        
        iso_5 = iso_4[2]  # get the isothermdata sliced at step 2
        @test num_steps(iso_5) == 1
        @test num_components(iso_5) == 3
        @test size(partial_pressures(iso_4)) == (7, 3)
        @test size(partial_pressures(iso_5)) == (1, 3)
        @test partial_pressures(iso_4)[2, :][2] == partial_pressures(iso_5)[2]
        @test partial_pressures(iso_4)[2, :] == partial_pressures(iso_5, step=1)

        @test partial_pressures(iso_4, step=2, component=:) == partial_pressures(iso_4)[2, :]
        @test partial_pressures(iso_5, component=1) == [2]
        
        iso_6 = iso_4[2, 1]
        @test partial_pressures(iso_4)[2, 1] == partial_pressures(iso_6)[1]

        iso_7 = IsothermData(partial_pressures_mpa = [[1, 2, 3], [4, 5, 6]])
        @test partial_pressures(iso_7, step=1) == [1, 4]
        
        iso_8 = IsothermData(partial_pressures_mpa = [[1, 2, 3]])
        @test partial_pressures(iso_8) == [1; 2; 3;;]
        iso_9 = IsothermData(partial_pressures_mpa = [1, 2, 3])
        iso_10 = IsothermData(partial_pressures_mpa = [1 2 3])
        @test partial_pressures(iso_8) == partial_pressures(iso_9) != partial_pressures(iso_10)
        
        # slicing
        @test num_steps(IsothermData()[0, 0]) == 0
        @test num_components(IsothermData()[0, 0]) == 0

        iso_11 = iso_4[1:3]
        @test partial_pressures(iso_11) == partial_pressures(iso_4, step=1:3)
        @test partial_pressures(iso_4, step=1:3, component=2:3) == partial_pressures(iso_4[1:3, 2:3])
        iso_12 = iso_4[3:4, 2:3]
        @test partial_pressures(iso_12, step=2) == partial_pressures(iso_4, step=4)[2:3]
        
        iso_13 = iso_4[:, :]
        iso_14 = iso_4[:]
        @test partial_pressures(iso_13) == partial_pressures(iso_14) == partial_pressures(iso_4)

        
        # removing desorption steps
        pres = [[1, 2, 3, 4, 5, 6, 7], [2, 3, 4, 5, 6, 7, 8]]
        conc = [[1, 2, 3, 4, 3, 2, 1], [10, 20, 30, 40, 30, 20, 30]]
        conc_invalid = [[1, 2, 3, 4, 3, 2, 1], [10, 20, 30, 40, 50, 20, 30]]
        iso_desorb = IsothermData(partial_pressures_mpa=pres, concentrations_cc=conc)
        iso_desorb_invalid = IsothermData(partial_pressures_mpa=pres, concentrations_cc=conc_invalid)
        @test ismissing(MembraneBase.increasing_concentration(iso_desorb_invalid; silent=true))
        
        # can we filter for only sorbing parts of the isotherm?
        iso_sorb = MembraneBase.increasing_concentration(iso_desorb)
        iso_sorb_redundant = MembraneBase.increasing_concentration(iso_sorb)
        @test partial_pressures(iso_sorb, component=1, step=1) == partial_pressures(iso_sorb_redundant, component=1, step=1)
        @test num_steps(iso_sorb) == 4 && num_components(iso_sorb) == 2
        @test partial_pressures(iso_sorb, step = 3) == [3, 4]
        allocs = @ballocated MembraneBase.increasing_concentration($iso_desorb) 
        @show allocs == 432
            
        # can we get the non-sorption components of the isotherm?
        iso_desorb_only = MembraneBase.remove_increasing_concentration_steps(iso_desorb)
        @test num_steps(iso_desorb_only) == 4 && num_components(iso_desorb_only) == 2   # coincidence
        @test partial_pressures(iso_desorb_only, step = 3) == [6, 7]
        allocs = @ballocated MembraneBase.remove_increasing_concentration_steps($iso_desorb) 
        @show allocs == 432

        # what if there's no valid point? We should get a single step
        iso_null = MembraneBase.increasing_concentration(iso_desorb_only)
        @test num_steps(iso_null) == 1
        @test num_components(iso_null) == 2

        # there was an example where a single desorption step didn't work
        pres = [[1, 2, 3, 4, 5]]
        conc = [[1, 2, 3, 4, 3]]
        iso_desorb_only_one = IsothermData(partial_pressures_mpa=pres, concentrations_cc=conc)
        iso_sorb = increasing_concentration(iso_desorb_only_one)
        @test num_steps(iso_desorb_only_one) == 5
        @test num_steps(iso_sorb) == 4
        iso_desorb_only = remove_increasing_concentration_steps(iso_desorb_only_one)
        @test num_steps(iso_desorb_only) == 2

        # BenchmarkTools allocations
        # todo optimize mole_fractions
        # allocs = @allocated mole_fractions(iso_3)
        # @show allocs
        # @show @btime mole_fractions($iso_3)
    end

    @testset "Isotherm Datasets - Deprecated" begin
        iso = IsothermData(partial_pressures_mpa = [1, 2, 3], concentrations_cc = [1, 4, 8], temperature_k = 273.15)
        iso2 = IsothermData(partial_pressures_mpa = [0.34, 2, 3, 8], concentrations_cc = [1, 2, 6, 7], temperature_k = 273.15)
        iso3 = IsothermData(partial_pressures_mpa = [0.4, 2, 5], concentrations_cc = [0.9, 4, 5], temperature_k = 300.15)
       
    end

    @testset "Transient Step Structure" begin
        # resampling
        tstep = TransientStepData([1., 2,   3,    4,     5,      6,    7,    8,    9,     10], 
                                [0, 0.5, 0.75, 0.875, 0.9375, 0.97, 0.98, 0.99, 0.995, 1.00])
        @test resample(tstep, 10, :Linear).time[6] == 6.0
        @test resample(tstep, 10, :Root).time[7] == 5.961012293408168
        @test resample(tstep, 5, :Log).time[3] == 3.1622776601683795

        tstep_err = TransientStepData([1 ± 3, 2,   3,    4,     5,      6,    7,    8,    9,     10], 
                                [0 ± 0.1, 0.5, 0.75, 0.875, 0.9375, 0.97, 0.98, 0.99, 0.995, 1.00])
        tstep_clean = strip_measurement_to_value(tstep_err)


        tstep_invalid = TransientStepData(
            [1 ± 0.1, 2,   1.9,    4,     5,      6,    7,    8,    9,     10], 
            [0 ± 0.2, 0.5, 0.75,   0.875, 0.9375, 0.97, 0.98, 0.99, 0.995, 1.00])
        @test_throws DomainError resample(tstep_invalid, 4, :Root)
        # @btime TransientStepData($[1, 2,   3 ± 0.1, 4, 5 ± 0.1, 6, 7, 8,  9, 10], 
        #     $[0, 0.5, 0.75, 0.875, 0.9375, 0.97, 0.98, 0.99, 0.995 ± 0.1, 1.00])
        
        ds = dataset(tstep)
        @test ds[2][2] == 0.5

        @test eltype(ds) <: Tuple{Float64, Float64}
        recovered_tsd = TransientStepData(ds)
        @test recovered_tsd.dimensionlesssorption[6] == 0.97
        # @btime resample($tstep, 10, :Root)
        # @btime dataset($tstep)
        
    end
end # end overall tests
return nothing