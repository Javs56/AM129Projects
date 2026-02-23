# file: complare_solutions.py
# Author: Sean Riedel

import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__=="__main__":
    # Get file name from command line
    if len(sys.argv)<2:
        print('No data file specified. Please specify one. For example:')
        print('python compare_solutions.py ./data/sol_021.dat')
        quit()
    fname = sys.argv[1]
    N = str(int(fname.split('_')[1].split('.')[0]))

    # Read in data file
    data = np.loadtxt(fname)

    # create exact solution data
    x = np.linspace(0, 2*np.pi, 500)
    exact_solution = np.exp(np.sin(x)) - np.exp(np.cos(3*x))

    # Create figure and plot
    fig,ax = plt.subplots()
    
    ax.plot(x,exact_solution)
    ax.scatter(data[:,0],data[:,1], color='red', marker='+')
    ax.legend(['Exact solution using ','Numerical solution'])
    ax.set_xlabel('x',fontsize='x-large')
    ax.set_ylabel('u(x)',fontsize='x-large')
    ax.set_title('Comparison of exact solution and numerical solution\nwith grid resolution of N='+N,fontsize='x-large')

    # Show figure
    plt.savefig("comparison_N"+N+".pdf")
    plt.show()
