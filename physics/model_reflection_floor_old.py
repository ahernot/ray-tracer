import numpy as np
from scipy import interpolate
import physics.model_velocity as VM
from physics.profile_temperature import calc_T
from physics.profile_pressure import calc_p


SPEED = {'simp' : VM.sound_velocity_simplified, 'leroy' : VM.sound_velocity_leroy, 'medwin' : VM.sound_velocity_medwin, 'mackenzie' : VM.sound_velocity_mackenzie}





def rho_seawater(S, z):
    """
    from : IES 80 – High Pressure International Equation of State of Seawater : https://unesdoc.unesco.org/ark:/48223/pf0000047363
    
    input :

    S: float : Practical salinity (‰), range:   0‰ ≤ S ≤ 42‰
    z: float : depth (m), range: -1000 ≤ z ≤ 0

    output

    rho : float : Seawater density (in kg.m^-3)
    """

    T = calc_T(z)

    #pressure linearily increases with depth
    p = calc_p(z)  # pressure in bar
 
    # Density of the Standard Mean Ocean Water (SMOW) [Bigg, 1967]
    rho_W = 999.842594 + 6.793952e-2 * T - 9.09529e-3 * np.power(T, 2) + 1.001685e-4 * np.power(T, 3) - 1.120083e-6 * np.power(T, 4) + 6.536336e-9 * np.power(T, 5)
    
    # Coefficients
    A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * np.power(T, 2) - 8.2467e-7 * np.power(T, 3) + 5.3875e-9 * np.power(T, 4)
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * np.power(T, 2)
    C = 4.8314e-4

    # One Atmosphere International Equation of State of Seawater (1980) (standard error = 3.6e-3 kg.m-3; validity 0‰≤S≤42‰, -2℃≤T≤40℃)
    rho_1atm = rho_W + A * S + B * np.power(S, 1.5) + C * np.power(S, 2)

    # Pure water terms
    K_W = 19652.21 + 148.4206 * T - 2.327105 * np.power(T, 2) + 1.360477e-2 * np.power(T, 3) - 5.155288e-5 * np.power(T, 4)
    A_W = 3.239908 + 1.43713e-3 * T + 1.16092e-4 * np.power(T, 2) - 5.79905e-7 * np.power(T, 3)
    B_W = 8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * np.power(T, 2)

    K_1atm = K_W + \
        (54.6746 - 0.603459 * T + 1.09987e-2 * np.power(T, 2) - 6.1670e-5 * np.power(T, 3)) * S + \
        (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * np.power(T, 2)) * np.power(S, 1.5)
    A = A_W + \
        (2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * np.power(T, 2)) * S + \
        1.91075e-3 * np.power(S, 1.5)
    B = B_W + \
        (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * np.power(T, 2)) * S

    # Secant Bulk Modulus
    K = K_1atm + A * p + B * np.power(p, 2)

    # High Pressure International Equation of State of Seawater
    rho = rho_1atm / (1 - p / K)
    
    return rho


def characteristic_specific_acoustic_impedance(density, speed) :
    Z = density*speed
    return Z

depth = np.arange(-1000,1,50)
rho = rho_seawater(35,depth)
mean_rho = np.mean(rho)
mean_c =1500

ZP0_REFERENCE_VALUES = {
    "estimation" : mean_rho*mean_c
} #unit : kg/(s*cm**2)


def compute_estimation_Zp0 ( speed_type = "simp" ) :

    assert speed_type in SPEED.keys(), f"speed_type must be among {set([key for key in SPEED.keys()])} and can't be \"{speed_type}\""

    depth = np.arange(-1000,1,50)
    c = SPEED[speed_type](35, calc_T(depth), depth)

    for indice in range(depth.shape[0]) :
        ZP0_REFERENCE_VALUES[depth[indice]] = c[indice]


#values from https://www.researchgate.net/publication/321661766_Acoustic_impedance_properties_of_seafloor_sediments_off_the_coast_of_Southeastern_Hainan
ZP1_REFERENCE_VALUES = {

    "silty clay" : 3160000,
    "sandy silt" : 2680000,
    "sand silt clay" : 2660000,
    "silt" : 2550000,
    "clayey silt" : 2390000
} 

estimation = np.mean([value for value in ZP1_REFERENCE_VALUES.values()])

ZP1_REFERENCE_VALUES["estimation"] = estimation
#unit : kg/(s*m**2)


ZP2_REFERENCE_VALUES = {
    "estimation" : 17000000
} 
#unit : kg/(s*m**2)

def  reflection_coefficient(theta0, Zp0 = ZP0_REFERENCE_VALUES["estimation"], Zp1 = ZP1_REFERENCE_VALUES["estimation"], Zp2 = ZP2_REFERENCE_VALUES["estimation"], D=2, Zs2 = None, thetas2 = None, absolute = False ) :
    """
    from https://cdn.intechopen.com/pdfs/45578/InTech-Ray_trace_modeling_of_underwater_sound_propagation.pdf (eq 21)

    input :

    theta0 : float : oriented angle between ray and the normal to bottom surface (grazing angle)
    D : float : sediment layer thickness (m), default is 2
    Zp0 : float : specific acoustic impedance for compression in the water column (kg/(s*m**2)), default is a mean value
    Zp1 : float : specific acoustic impedance for compression in the sediment layer (kg/(s*m**2)), default is a mean value
    Zp2 : float : specific acoustic impedance for compression in the solid half-space (kg/(s*m**2)), default is a mean value
    Zs2 : float : specific acoustic impedance for shear in the solid half-space (kg/(s*m**2)), default is None ==> shear is not considered
    thetas2 : float : grazing angle of the transmitted shear wave in the solid half-space, default is None ==> shear is not considered
    absolute : bool : if True, the function will return the absolute value of the coefficient
    
    output :

    if absolute :
        Rb : numpy.complex128 : reflection coefficient
    if not absolute :
        rb : numpy.float64 : absolute value of the reflection coefficient
    """
        
    gammap1 = np.sin(theta0)

    a = 2 * 1j * gammap1 * D
    
    r01 = (Zp1-Zp0)/(Zp1+Zp0)
    if Zp2 is None or thetas2 is None : #shear is not considered
        r12 = (Zp2-Zp1)/(Zp2+Zp1)
    else : #if shear is considered
        r12 = (Zp2*np.cos(thetas2)**2+Zs2*np.sin(thetas2)**2-Zp1)/(Zp2*np.cos(thetas2)**2+Zs2*np.sin(thetas2)**2+Zp1)
    
    Rb = (r01+r12*np.exp(-a))/(1+r01*r12*np.exp(-a))
    
    if absolute :
        rb = np.abs(Rb)
        return rb
    else :
        return Rb

