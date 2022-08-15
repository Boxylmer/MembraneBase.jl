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

        # the inverse hessian methods
        objective_function(xy) = ((xy[1] + xy[2])^2 + (xy[1]-1)^2 + (xy[2]-0.8)^2)^2
        minimizer = [0.494234234, 0.25034623146]
        @test errs = rss_minimizer_standard_errors(objective_function, minimizer, 10) ≈ [0.15049509481360718, 0.1517352228363271]

    end

    @testset "Statistical Methods" begin
        ndata = 50
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
        @test round(j_σ[1]; digits=0) == 1

        b_σ = bootstrap_uncertainty(fitter, data)
        @test round(b_σ[1]; digits=0) == 1

        function run_boot()
            b_σ = bootstrap_uncertainty(fitter, data)
        end
        function run_jack()
            j_σ = jackknife_uncertainty(fitter, data)
        end
        # @btime $run_boot()
        # @btime $run_jack()
        

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
        @test mole_fractions(iso_2) == [1; 1; 1][:,:]
        
        # BenchmarkTools allocations
        allocs = @allocated mole_fractions(iso_3)
        @show allocs
        @show @btime mole_fractions($iso_3)


    end

    @testset "Isotherm Datasets - Deprecated" begin
        iso = IsothermData(partial_pressures_mpa = [1, 2, 3], concentrations_cc = [1, 4, 8], temperature_k = 273.15)
        iso2 = IsothermData(partial_pressures_mpa = [0.34, 2, 3, 8], concentrations_cc = [1, 2, 6, 7], temperature_k = 273.15)
        iso3 = IsothermData(partial_pressures_mpa = [0.4, 2, 5], concentrations_cc = [0.9, 4, 5], temperature_k = 300.15)
        # single_isotherm_dataset = MembraneBase.TPCDataset(iso)
        # @test single_isotherm_dataset[2] == [273.15, 2.0, 4.0]
        # @test single_isotherm_dataset[1:2] == [[273.15, 1.0, 1.0], [273.15, 2.0, 4.0]]
        # new_data = copy(single_isotherm_dataset)
        # @test deleteat!(new_data, 2) == [[273.15, 1.0, 1.0], [273.15, 3.0, 8.0]]

        # multi_dataset = MembraneBase.TPCDataset([iso, iso2, iso3])
        # @test multi_dataset.tpc_vectors == [
        #     [273.15, 0.34, 1.0],
        #     [273.15, 1.0, 1.0],
        #     [273.15, 2.0, 2.0],
        #     [273.15, 2.0, 4.0],
        #     [273.15, 3.0, 6.0],
        #     [273.15, 3.0, 8.0],
        #     [273.15, 8.0, 7.0],
        #     [300.15, 0.4, 0.9],
        #     [300.15, 2.0, 4.0],
        #     [300.15, 5.0, 5.0]]
        
        # recovered_isotherms = MembraneBase.get_isotherms(multi_dataset)
        
        # remade_dataset = MembraneBase.TPCDataset(recovered_isotherms)
        # @test remade_dataset == multi_dataset

        # test the potential edge case of one last temperature remaining that's different
        # iso4 = IsothermData(partial_pressures_mpa = [1, 2, 3], concentrations_cc = [1, 4, 8], temperature_k = 273.15)
        # iso5 = IsothermData(partial_pressures_mpa = [6], concentrations_cc = [6], temperature_k = 308.15)
        # iso6 = IsothermData(partial_pressures_mpa = [5], concentrations_cc = [5], temperature_k = 100.15)
        # one_single_temp_dataset = MembraneBase.TPCDataset([iso4, iso5, iso6])
        # recovered_isotherms = MembraneBase.get_isotherms(one_single_temp_dataset)
        # remade_dataset = MembraneBase.TPCDataset(recovered_isotherms)
        # @test remade_dataset == one_single_temp_dataset

        # test the potential edge case of a single datapoint
        # iso7 = IsothermData(partial_pressures_mpa = [1], concentrations_cc = [1], temperature_k = 273.15)
        # dataset = MembraneBase.TPCDataset([iso7])
        # recovered_isotherms = MembraneBase.get_isotherms(dataset)
        # remade_dataset = MembraneBase.TPCDataset(recovered_isotherms)
        # @test remade_dataset == dataset

        # toss all the isotherms together for good measure
        # dataset = MembraneBase.TPCDataset([iso, iso2, iso3, iso4, iso5, iso6, iso7])
        # recovered_isotherms = MembraneBase.get_isotherms(dataset)
        # remade_dataset = MembraneBase.TPCDataset(recovered_isotherms)
        # @test remade_dataset == dataset


        # todo the same thing as above, but with use_fugacity set true

        # fug_iso = IsothermData(fugacities_mpa = [1, 2, 3], concentrations_cc = [1, 4, 8], temperature_k = 273.15)
        # fug_dataset = MembraneBase.TPCDataset(fug_iso; use_fugacity = true)
        # @test fugacities(IsothermData(dataset); component=1, step = 3) == 3
        # @show IsothermData(dataset)
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
nothing