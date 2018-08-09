#==========================================================================================================#
# DENSITY_CALC.PY
#==========================================================================================================#

#!/opt/local/bin/python

# To calculate densities given a mass and a radius

#==========================================================================================================#

# Import modules

import numpy as np
import math
# Set global variables

mu         = 2.8 # Mean mass per particle
mh         = 1.6737236e-27 # Mass of a hydrogen atom in kg
msun       = 1.989e30 # Mass of sun in kg
pc         = 3.08567758e16 # A parsec in metres
percm2perm = 1.e6

def column_density(Flux, Beam, Wave, Temp, Kappa, gas2dust):

    from planck_func import planck_wave

    B = planck_wave( Wave, Temp )

    Omega = (math.pi / (4.0 * math.log(2.)) ) *Beam[0]*Beam[1]

    N = (Flux * gas2dust) / (Omega * mu * (mh*1.e3) * Kappa * B)

    return N

# Define subroutines

def number_density_sphere_pc( Mass_sol, Radius_pc, mu ):

    # This subroutine accepts mass in solar masses and radius in pc and calculates the number density.

    Mass = Mass_sol * msun
    Radius = Radius_pc * pc

    n = Mass / (((4. / 3.)*np.pi) * mu * mh * Radius**3.0)

    # Convert to astronomical units particles per cubic centimetre

    n = n / percm2perm

    return n

def mass_density_sphere( Mass_sol, Radius_pc ):

    # This subroutine accepts mass in solar masses and radius in pc and calculates the mass density.

    Mass = Mass_sol * msun
    Radius = Radius_pc * pc

    rho = Mass / (((4. / 3.)*np.pi) * Radius**3.0)

    return rho
