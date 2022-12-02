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
const PA_PER_MPA = 1 / MPA_PER_PA 
const MPA_PER_BAR = 0.1
const PA_PER_ATM = MPA_PER_ATM / MPA_PER_PA
const ATM_PER_PA = 1 / PA_PER_ATM
const MPA_PER_PSI = 0.006894757293
const PSI_PER_MPA = 1 / MPA_PER_PSI
const ATM_PER_PSI = ATM_PER_MPA * MPA_PER_PSI
const PSI_PER_ATM = 1 / ATM_PER_PSI

# volume
const L_PER_CC = 0.001

# moles
const CC_PER_MOL_STP = 22414

# permeability
const BARRER_PER_CC_CC_MPA_CM2_S = 1 / (750.06156130264 * 1e-10)
const CC_CC_MPA_CM2_S_PER_BARRER = 1 / BARRER_PER_CC_CC_MPA_CM2_S
const BARRER_PER_CCSTP_CM_S_CMHG = 1e-10   # Barrer / (CCstp / (cm2 * s * cmHg))
const CCSTP_CM_S_CMHG_PER_BARRER = 1 / BARRER_PER_CCSTP_CM_S_CMHG
const BARRER_PER_CCSTP_CM_S_ATM = 75.99998769899 * BARRER_PER_CCSTP_CM_S_CMHG  # Barrer / (CCstp / (cm2 * s * atm))
const CCSTP_CM_S_ATM_PER_BARRER = 1 / BARRER_PER_CCSTP_CM_S_ATM
const BARRER_PER_MOL_CM_S_ATM = BARRER_PER_CCSTP_CM_S_ATM / CC_PER_MOL_STP  # Barrer / (mol / (cm2 * s * atm))
const MOL_CM_S_ATM_PER_BARRER = 1 / BARRER_PER_MOL_CM_S_ATM