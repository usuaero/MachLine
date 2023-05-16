import flow54 as fl
import numpy as np

study_dir = "studies/supersonic_double_wedge_wing/"


def get_se_pressures(M, alpha, half_angle):
    # Runs the comparison of the diamond wing to shock-expansion theory

    # Parameters
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    rho = 1.225
    p_inf = 1.0e5
    
    # Run shock-expansion comparison
    airfoil = fl.DiamondAirfoil(half_angle, 1.0)
    airfoil.set_state(M, alpha, gamma, p_inf, T_inf, c_inf, rho, 1.0e-5, 800.0, 0.7)
    p2, p3, p4, p5 = airfoil.get_pressures()
    x = 0.5*gamma*M**2
    Cp2 = (p2/p_inf-1.0)/x
    Cp3 = (p3/p_inf-1.0)/x
    Cp4 = (p4/p_inf-1.0)/x
    Cp5 = (p5/p_inf-1.0)/x

    return Cp2, Cp3, Cp4, Cp5


if __name__=="__main__":

    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)
    half_angle = 5

    Cps = np.zeros((len(Ms), len(alphas), 4))

    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):
            Cps[j,k,:] = get_se_pressures(M, alpha, half_angle)

    np.savetxt(study_dir + "shock_expansion_data.csv", Cps.reshape((len(Ms)*len(alphas), 4)), delimiter=',')