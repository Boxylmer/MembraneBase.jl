
"""
    IsothermData(; 
        partial_pressures_mpa=nothing, concentrations_cc=nothing, activities=nothing, temperature_k=nothing, 
        rho_pol_g_cm3=nothing, pen_mws_g_mol=nothing,
    )
Container for all data relevant to isotherms.

Data specific to each component and step are organized in matrix format with columns representing components and rows representing steps
i.e., with s -> step and c -> component:
```
    [ s1c1  s1c2  s1c3 ...
      s2c1  s2c2  s2c3 ...
      s3c1  s3c2  s3c3 ...
      ...   ...   ...  ... ]
```

# Inputs 
Inputs can be vectors of vectors e.g., `partial_pressures_mpa = [[component_1-step_1, component_1-step_2, ...], [component_2-step_1, component_2-step_2, ...], ...]`
Or, if only one component is present, only a vector e.g., `[step_1, step_2, ...]`

"""
struct IsothermData{PPT, CT, AT, FT, TT, PDT, PMWT}
    partial_pressures::PPT              # MPa
    concentrations::CT                  # CC(STP)/CC(Pol)
    activities::AT                      # Unitless
    fugacity::FT                        # MPa
    temperature::TT                     # K
    polymer_density::PDT                # g/cm3
    penetrant_molecular_weights::PMWT   # g/mol
    num_components::Int64
    num_steps::Int64
end

module IsothermHelperFunctions
    """Convert a vector of vectors (vector of components holding vectors of steps) into a [step, component] matrix. 
        Matrices are passed through and vectors are assumed to be one component"""
    function convert_data_input_to_formatted_matrix(vec_or_vec_of_vecs_or_matrix)
        if !isnothing(vec_or_vec_of_vecs_or_matrix)
            # pesky tuples cause hcat to mess up
            if typeof(vec_or_vec_of_vecs_or_matrix) <: Tuple
                return convert_data_input_to_formatted_matrix([item for item in vec_or_vec_of_vecs_or_matrix])
            end

            # handle element types meanings
            if eltype(vec_or_vec_of_vecs_or_matrix) <: Number
                vec_of_vecs = [vec_or_vec_of_vecs_or_matrix] 
            elseif eltype(vec_or_vec_of_vecs_or_matrix) <: Matrix
                return vec_or_vec_of_vecs_or_matrix
            elseif eltype(vec_or_vec_of_vecs_or_matrix) <: Vector
                vec_of_vecs = vec_or_vec_of_vecs_or_matrix
            
            else
                throw(ArgumentError("Review the datatypes passed to convert_vec_to_matrix"))
            end
            resultant_matrix = hcat(vec_of_vecs...)
        else
            resultant_matrix = nothing
        end
    end

    function convert_data_input_to_formatted_vector(num_or_vec)
        if !isnothing(num_or_vec)
            # handle element types meanings
            if typeof(num_or_vec) <: Number
                return [num_or_vec]
            elseif eltype(num_or_vec) <: Number
                return num_or_vec 
            else
                throw(ArgumentError("Review the datatypes passed to this function"))
            end
        else
            return nothing
        end
    end

    """
    Get the step and component of some data matrix in the isotherm. 
    - Specifying neither `step` nor `component` returns the whole matrix.
    - Specifying one of the above returns a vector.
    - Specifying both will return a single value.
    """
    function get_data(data::AbstractArray; step = : , component = : )
        return @view data[step, component]
    end

    function get_data(data::Nothing; kwargs...)
        return nothing
    end

end  # end IsothermHelperFunctions


function IsothermData(; 
        partial_pressures_mpa=nothing, concentrations_cc=nothing, activities=nothing, fugacities_mpa=nothing, 
        temperature_k=nothing, rho_pol_g_cm3=nothing, pen_mws_g_mol=nothing,
    )
    # partial_pressures = [[component_1-step_1, component_1-step_2, ...], [component_2-step_1, component_2-step_2, ...], ...]
    # same formatting for concentrations. If you feed in [step_1, step_2, ...] for either, it's assumed there's only one component
    # temperature only needs to be a number
    partial_pressure_matrix = IsothermHelperFunctions.convert_data_input_to_formatted_matrix(partial_pressures_mpa)
    concentration_matrix = IsothermHelperFunctions.convert_data_input_to_formatted_matrix(concentrations_cc)
    activity_matrix = IsothermHelperFunctions.convert_data_input_to_formatted_matrix(activities)
    fugacity_matrix = IsothermHelperFunctions.convert_data_input_to_formatted_matrix(fugacities_mpa)
    molar_mass_vector = IsothermHelperFunctions.convert_data_input_to_formatted_vector(pen_mws_g_mol)

    num_steps, num_components = ensure_matrices_are_same_size(partial_pressure_matrix, concentration_matrix, activity_matrix)
    
    IsothermData(
        partial_pressure_matrix, concentration_matrix, activity_matrix, fugacity_matrix,
        temperature_k, rho_pol_g_cm3, molar_mass_vector, num_components, num_steps
    )
end

function Base.getindex(iso::IsothermData, step, component=:)
    nsteps = typeof(step) <: Colon ? num_steps(iso) : length(step) 
    ncomps = typeof(component) <: Colon ? num_components(iso) : length(component)
    return IsothermData(
        materialize(partial_pressures(iso; component, step)),
        concentration(iso; component, step),
        activities(iso; component, step),
        fugacities(iso; component, step),
        temperature(iso),
        polymer_density(iso),
        molecular_weights(iso, component),
        ncomps,
        nsteps,
    )
end

"""
    isotherm_dataset(x::Matrix, y::Matrix, component::Int)
converts two isotherm components into a vector steps in the format of [(type_x, type_y)...]
For multicomponent isotherms, the format is dataset[step] = [[xcomp1, xcomp2, ...], [ycomp1, ycomp2, ...]]
    e.g., `mydataset = Isotherm.dataset(myisotherm.partial_pressures, myisotherm.concentrations)` gives a vector of vectors in the format of:
        myset[1] = [[(pressure of component 1), (pressure of component 2), ...], [(conc of component 1), (conc of component 2), ...]]
"""
function isotherm_dataset(x::AbstractMatrix, y::AbstractMatrix, component::Int)
    data_x = x[:, component]
    data_y = y[:, component]
    return [collect(item) for item in zip(data_x, data_y)]
end 

function isotherm_dataset(x::AbstractMatrix, y::AbstractMatrix)
    return [[x[idx, :], y[idx, :]] for idx in 1:size(x)[1]]
end 

function strip_measurement_to_value(iso::IsothermData)  # strip all measurements from the isotherm and return a copy
    return IsothermData(
        strip_measurement_to_value(iso.partial_pressures), 
        strip_measurement_to_value(iso.concentrations), 
        strip_measurement_to_value(iso.activities), 
        strip_measurement_to_value(iso.fugacity),
        strip_measurement_to_value(iso.temperature), 
        strip_measurement_to_value(iso.polymer_density),
        strip_measurement_to_value(iso.penetrant_molecular_weights), 
        iso.num_components, 
        iso.num_steps,
    )
end

# define a whole bunch of nice getter functions
"""
    pressure(isotherm::IsothermData; step)
Get the total pressure at a step in the isotherm. If no step is provided, a vector of steps (`pressure[step]`) is returned.
"""
function pressure(isotherm::IsothermData; step=:)
    if typeof(step) <: Colon
        return sum(IsothermHelperFunctions.get_data(isotherm.partial_pressures; step=step); dims=2)[:]
    elseif typeof(step) <: Integer
        return sum(IsothermHelperFunctions.get_data(isotherm.partial_pressures; step=step))
    end
end
"""
    partial_pressures(isotherm::IsothermData; component, step)
Get the partial pressures of some component in the isotherm. Also synonymous with `pressures`.
"""
partial_pressures(isotherm::IsothermData; component=:, step=:) = IsothermHelperFunctions.get_data(isotherm.partial_pressures; component=component, step=step)
pressures = partial_pressures

"""
    concentration(isotherm::IsothermData; component=nothing, gas_units=:cc, pol_units=:cc)
Get the concentrations of a component in the isotherm with specific units. 
# Arguments
- `gas_units`: Units of the penetrant.
    - Can specify: `:g`, `:cc`
!!! note
    `gas_units` should really be `penetrant_units`, since liquid phase isotherms exist. For now, consider them synonyms.
- `polymer_units`: Units of the polymer.
    - Can specify: `:g`, `:cc`

For example `concentration(some_isotherm; gas_units=:g, pol_units=:cc)` will get the component concentration in units of ``\\frac{g}{CC_{polymer}}``

!!! warning
    This method is far from optimized. In fact, every time `concentration` is called, all isotherm components are converted to whatever units were specified. So right now, if performance is critical and all concentrations are required, just get them all at once and iterate as needed.

"""
function concentration(isotherm::IsothermData; component=:, step=:, gas_units=:cc, pol_units=:cc)

    # units start out in cc(STP)/cc(pol)
    if gas_units == :cc && pol_units == :cc
        return IsothermHelperFunctions.get_data(isotherm.concentrations; component=component, step=step)
    else
        returned_data = copy(isotherm.concentrations)
        
        if gas_units == :g
            if isnothing(molecular_weights(isotherm))
                throw(MissingException("Molecular weights were not provided."))
            end
            returned_data = returned_data ./ CC_PER_MOL_STP
            returned_data = returned_data .* transpose(molecular_weights(isotherm))
        end

        if pol_units == :g
            if isnothing(polymer_density(isotherm))
                throw(MissingException("Polymer density was not provided."))
            end
            returned_data = returned_data ./ polymer_density(isotherm)
        end

        return IsothermHelperFunctions.get_data(returned_data; component=component, step=step)
    end
end

"""
    mole_fractions(isotherm::IsothermData, step=:)
Get the mole fractions of some component in the isotherm.
Mole fractions are assumed to be partial pressures over total pressure.  
"""
function mole_fractions(isotherm::IsothermData; step=:)
    pressures = pressure(isotherm; step=step)
    partial_pressures = IsothermHelperFunctions.get_data(isotherm.partial_pressures; step=step)
    return partial_pressures ./ pressures
end

"""
    activities(isotherm::IsothermData, component=:, step=:)
Get the activities of some component in the isotherm. 
"""
activities(isotherm::IsothermData; component=:, step=:) = IsothermHelperFunctions.get_data(isotherm.activities; component=component, step=step)

"""
    fugacity(isotherm::IsothermData, component=:, step=:)
Get the fugacity of some component in the isotherm. 
"""
fugacities(isotherm::IsothermData; component=:, step=:) = IsothermHelperFunctions.get_data(isotherm.fugacity; component=component, step=step)


"""
    temperature(isotherm::IsothermData)
Get the isotherm temperature. 
"""
temperature(isotherm::IsothermData) = isotherm.temperature


"""
    polymer_density(isotherm::IsothermData)
Get the density of the polymer in the isotherm. 
"""
polymer_density(isotherm::IsothermData) = isotherm.polymer_density


"""
    molecular_weights(isotherm::IsothermData, component=nothing)
Get molecular weights of a component in the isotherm. 
"""
function molecular_weights(isotherm::IsothermData, component=nothing)
    if isnothing(isotherm.penetrant_molecular_weights)
        return nothing
    elseif !isnothing(component)
        return isotherm.penetrant_molecular_weights[component]
    else
        return isotherm.penetrant_molecular_weights
    end
end


"""
    num_components(isotherm::IsothermData)
Get the number of components present in the isotherm. 
"""
num_components(isotherm::IsothermData) = isotherm.num_components


"""
    num_steps(isotherm::IsothermData)
Get the number of steps present in the isotherm. 
"""
num_steps(isotherm::IsothermData) = isotherm.num_steps


# getter functions that do downstream caclulations on the isotherm data

"""
    mass_sorbed(isotherm::IsothermData, polymer_mass_g::Number, component=:)
Get mass (in **g**) of penetrant sorbed into the polymer phase, given the mass (in **g**) of the polymer.
!!! note
    If no component is specified in a multicomponent isotherm, a `vector` of masses is returned, **not** the total mass sorbed.  
!!! warning
    *The mass of the polymer* in this case means the mass of the polymer *only*. Not the mass of the polymer phase.
    In most cases, this distinction is insignificant, however in very highly sorbing polymers, this difference can cause significant error if not considered. 
"""
function mass_sorbed(isotherm::IsothermData, polymer_mass_g::Number; component=:, step=:)
    concs_g_g = concentration(isotherm; component=component, step=step, gas_units=:g, pol_units=:g)
    _masses_sorbed = concs_g_g .* polymer_mass_g
    return _masses_sorbed
end

# NOT FULLY GENERIC YET: assumes that all other components are ≈ 0. Could circumvent this by adding the masses of the other components into the denominator
function penetrant_mass_fractions(isotherm::IsothermData; component=:, step=:)
    concs_g_g = concentration(isotherm; component=component, step=step, gas_units=:g, pol_units=:g)
    all_concs_in_step = concentration(isotherm;step=step, gas_units=:g, pol_units=:g)
    mass_fractions = concs_g_g ./ (1 .+ sum(all_concs_in_step))  # system = 1g pol + whatever g penetrant

    return mass_fractions

end