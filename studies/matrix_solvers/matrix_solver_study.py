import os
import subprocess as sp
import numpy as np

if __name__=="__main__":

    # Solver options
    solvers = ["GMRES", "QRUP", "BJAC", "BSOR", "LU"]