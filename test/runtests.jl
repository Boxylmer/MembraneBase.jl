using MembraneUtils
using Test

@testset "MembraneUtils.jl" begin
    
    # ensure_matrices_are_same_size
    mat1 = rand(1, 2, 3)
    mat2 = rand(1, 2, 3)
    mat3 = rand(1, 3, 3)
    @test ensure_matrices_are_same_size(mat1, mat2) == size(mat1)
    @test_throws ensure_matrices_are_same_size(mat1, mat3)

    # fitting linear data
    x = [1, 2, 3]
    y = [2, 4.1, 5.8]
    y_meas = [2 ± 0.0, 4.1 ± 0.4, 5.8 ± 0.3]
    @test_logs (:warn,) PolymerMembranes.fit_linear_data(x, y_meas)  # using a measurement with perfect precision should throw a warning
    @test_nowarn PolymerMembranes.fit_linear_data(x, y)  

    # conversion methods
    molecular_weights = [1, 3., 4.13, π, 1e30]
    mole_fractions = [0.2, 0.4, 0.1, 0.2, 0.1]

    mass_fractions = PolymerMembranes.mole_fractions_to_mass_fractions(mole_fractions, molecular_weights)
    recovered_mole_fractions = PolymerMembranes.mass_fractions_to_mole_fractions(mass_fractions, molecular_weights)
    @test mole_fractions ≈ recovered_mole_fractions

    molar_volume = 0.228  # L/mol
    density = PolymerMembranes.molar_volume_to_density(molar_volume, mole_fractions, molecular_weights)
    recovered_molar_volume = PolymerMembranes.density_to_molar_volume(density, mole_fractions, molecular_weights)
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
    taylorfunc = PolymerMembranes.get_taylor_series_function(myfunc, [0.5, 0.5, 0.5, 0.5, 0.5], 17)
    @test taylorfunc([0.1, 0.5, 0.1, 0.2, 0.56]) ≈ myfunc([0.1, 0.5, 0.1, 0.2, 0.56])
   


end
