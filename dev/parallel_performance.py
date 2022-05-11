import os

import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt


if __name__=="__main__":

    # Declare desired numbers of threads
    Ns = [1, 2, 4, 8]
    N_test = 3 # Number of test runs for averaging purposes

    # Initialize time storage
    times = np.zeros((len(Ns), N_test))

    # Get current environment variables
    old_env = os.environ.copy()

    # Loop through numbers of threads
    for i, N_threads in enumerate(Ns):

        # Set number of threads
        os.environ.update({"OMP_NUM_THREADS" : str(N_threads)})

        # Loop through tests
        for j in range(N_test):

            # Run py.test
            result = sp.run(['py.test'], capture_output=True, text=True)

            # Check for fail
            if 'failed' in result.stdout:
                print(result.stdout)

            # Determine how long it took
            output_lines = result.stdout.splitlines()
            time = output_lines[-1].split()[-2]
            times[i,j] = float(time.replace('s', ''))

    # Reset environment variables
    os.environ.clear()
    os.environ.update(old_env)

    # Compute average times
    avg_times = np.sum(times, axis=1)/N_test

    # Plot
    plt.figure()
    plt.plot(Ns, avg_times, 'kx-')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of Threads')
    plt.ylabel('Avg. Computation Time [s]')
    plt.show()