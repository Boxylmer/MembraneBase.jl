const R_J_MOL_K = 8.31446261815324            # J / (K * mol)
const R_ATM_L_K_MOL = 0.0820573660809596      # atm * L / (K * mol)

const R_PA_M3_K_MOL = R_J_MOL_K               # Pa * m^3 / (K * mol) 
const R_MPA_L_K_MOL = R_J_MOL_K / 1000        # MPa * L / (K * mol)
const R_MPA_CM3_K_MOL = R_MPA_L_K_MOL * 1000  # MPa * cm^3 / (K * mol)
const R_BAR_M3_K_MOL = R_J_MOL_K / 100000     # Bar * m^3 / (K * mol)
const R_ATM_M3_K_MOL = R_ATM_L_K_MOL / 1000   # atm * m^3 / (K * mol)

# pressure
const MPA_PER_ATM = 0.101325
const ATM_PER_MPA = 1 / MPA_PER_ATM 
const MPA_PER_PA = 0.000001
const MPA_PER_BAR = 0.1
const PA_PER_ATM = MPA_PER_ATM / MPA_PER_PA
const ATM_PER_PA = 1 / PA_PER_ATM

# volume
const L_PER_CC = 0.001

const CC_PER_MOL_STP = 22414

# permeability
const BARRER_PER_CC_CC_MPA_CM2_S = 1 / (750.06156130264 * 1e-10)
const CC_CC_MPA_CM2_S_PER_BARRER = 1 / BARRER_PER_CC_CC_MPA_CM2_S