
# Convert Draine & Li model structures to SEDs with flux normalized
# to something.  Here I'm trying to normalize to 1 solar mass of dust
# at 10 Mpc.  This should work on either a single model dictionary
# or a list of models.

# DL models give: nu P_nu in erg/s/H nucleon
#     j_nu in Jy cm^2 /sr/H nucleon
# make the output SEDs in Jy at 10 Mpc per 1 Msun of dust
# 1 Jy = 1e-23 erg/s/cm^2/Hz
# if I use j_nu, I have Jy in and Jy out, so just need to scale
# by the distance in cm squared, and the mass of dust
# (H nucleons per Msun * dust/gas ratio)


# bjw 5/28/2016

# DL models have jnu per H nucleon
# M_dust = jnu * (M_H / m_H) * (M_gas / M_H) * (M_dust / M_gas)
#  factors are: Msun / m_H, Mgas/M_H = 1.36 for He, and
#  dust/H ~0.002-0.01 dep. on SMC/LMC/MW.

# flux jnu = M_dust * (m_H / M_H) * (M_H / M_dust)
# I think the m_H / M_H factor should really be n_H per Msun, which is
# the reciprocal; possible to get confused here but to get flux
# per Msun, since the original number is flux per proton,
# clearly you have to scale up by protons per solar mass
# and to get flux per Msun of dust, you must scale up by ~100x
# since there is the same flux but 100x less mass of dust.

# Dust to H ratios come from Table 3 of Draine & Li 2007
# models are 1-11
# name: MW3.1_10, MW3.1_20, MW3.1_30, MW3.1_40, MW3.1_50, MW3.1_60,
#       LMC2_00, LMC2_05, LMC2_10, SMC
# dust/H: 0.0100, 0.0100, 0.0101, 0.0102, 0.0102, 0.0103, 0.0104,
#         0.00343, 0.00344, 0.00359, 0.00206
# This is dust/H ,and dust/gas = 1/1.36 * dust/H


def convert_draineli_sed(dlmodel):

    dust_to_h_array = [0.0100, 0.0100, 0.0101, 0.0102, 0.0102, 0.0103, 0.0104,
                   0.00343, 0.00344, 0.00359, 0.00206]
    nmod = len(dlmodel)
    sedstruct = []
    # 10 Mpc in cm
    refdist = 10.0 * 3.086e24
    atoms_per_msun = 1.989e33 / 1.67e-24
    # dust_to_gas = 0.001
    gas_per_h = 1.36
    for i in range(nmod):
        label = dlmodel[i]['uavg']
        name = dlmodel[i]['filename']
        wave = dlmodel[i]['sedwave']
        flux = dlmodel[i]['sedjnu']
        imodel = int(dlmodel[i]['grainmodel'])
        if imodel <1 or imodel>11:
            dust_to_h = 0.01
        else:
            dust_to_h = dust_to_h_array[imodel-1]
        # If over/underflow becomes a problem, do the calculation in log.
        #flux = flux * atoms_per_msun * gas_to_h * dust_to_gas / refdist**2
        factor = atoms_per_msun / dust_to_h / refdist**2
        flux = flux * factor
        sed1 = {'label': label, 'name':name, 'wave':wave, 'flux':flux}
        sedstruct.append(sed1)

    return sedstruct
