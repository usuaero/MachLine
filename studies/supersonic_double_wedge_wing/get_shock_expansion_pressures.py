import flow54 as fl
import numpy as np
import matplotlib.pyplot as plt

study_dir = "studies/supersonic_double_wedge_wing/"


def run_comparison(M, alpha, half_angle):
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

    grids = ["coarse", "medium", "fine", "ultra_fine"]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)
    half_angles = [5]

    Cps = np.zeros((len(grids), len(Ms), len(alphas), len(half_angles), 4))

    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                for l, half_angle in enumerate(half_angles):

                    Cps[i,j,k,l] = run_comparison(M, alpha, half_angle)

    np.savetxt(study_dir + "shock_expansion_data.csv", Cps.reshape((len(grids)*len(Ms)*len(alphas)*len(half_angles), 4)), delimiter=',')