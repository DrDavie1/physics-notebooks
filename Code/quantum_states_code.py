


# This file is used to demonstrate the functions used in the notebook assignment.
#--------------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
import scipy.special as sps

# input x array, output: harmonic potential array

def potential_qho(x: np.ndarray) -> np.ndarray:

    v = x**2
    
    return v

# plot the harmonic potential

# thoughout the code plotting functions have no return type

def plot_qho_potential(x_min: float,x_max: float,N: int) -> None:
    
    x_values = np.linspace(x_min,x_max,N)
    
    v_values= potential_qho(x_values)

    fig, ax = plt.subplots()

    ax.set(xlabel='$x$',ylabel='$V(x)$',title='Harmonic Potential')

    # scatter plot used to demonstrate discrete nature of code

    ax.scatter(x_values,v_values,color='black',s=5)

    ax.plot(x_values,v_values,color='red')


# input: x array and potential array, output: H matrix

def hamiltonian(x: np.ndarray,v: np.ndarray) -> np.ndarray:

    #construct matricies

    N = len(x)

    V = np.zeros((N,N))
    D = np.zeros((N,N))

    dx = x[1] - x[0] # assumes equal spacing (which is fine)


    for i in range(N):

        #Defining the potential matrix from array V(x_i)

        V[i][i] = v[i]

        #Defining D matrix

        D[i][i] = 2/((dx)**2)

        if i-1 >= 0:
            D[i][i-1] = -1/((dx)**2)
            
        if i+1 < N:   
            D[i][i+1] = -1/((dx)**2)

    #add matricies and return H

    H = D + V
        
    return H


# input: x range, N points, output: eigen values, eigen vectors, x values (useful for plotting)

def harmonic_eigen_val_vec(x_min: float,x_max: float,N: int) -> tuple:
    
    x_values = np.linspace(x_min,x_max,N)
    
    v_values = potential_qho(x_values)

    # define hamiltonian
    
    H = hamiltonian(x_values,v_values)

    # compute eigen values and vectors using scipy function
    
    eigen_values, eigen_vectors = sp.linalg.eigh(H)

    return eigen_values, eigen_vectors, x_values


# input: Number of expected values to compute, output: expected values for harmonic oscillator

def expected_harmonic_eigen_val(N_values: int) -> np.ndarray:

    # defined as the first N_values odd numbers, consequence of explanation in markdown cell 

    expected = np.linspace(1,(2*N_values - 1),N_values)

    return expected


import pandas as pd # for table

# input: eigen values, expected eigen values, output: table to compare values including % error 

def eigen_vector_table(eigen_values,expected_values) -> pd.DataFrame:

    percentage_error = abs((1 - (eigen_values/expected_values)) * 100)

    data = {f'Calculated':eigen_values,f'Expected':expected_values,'Error (%)':percentage_error}

    table = pd.DataFrame(data)

    return table


def probability_density(i: int,eigen_vectors: np.ndarray) -> np.ndarray:

    # select ith column and find modulus squared. 

    psi_squared = np.conj(eigen_vectors[:,i]) * eigen_vectors[:,i] 

    return psi_squared 


# plot wave functions for a given range of n values (0 -> n_tot)

def plot_harmonic_wavefunction(n_tot: int, x_values: np.ndarray, eigen_vectors: np.ndarray) -> None:

    # define size of axis based on chosen n_tot
    
    fig, axs = plt.subplots(n_tot,1,figsize=(6,4*n_tot))

    fig.suptitle(f'Harmonic Wavefunctions, n = [0,{n_tot-1}]')

    for n in range(0,n_tot):

        prob = probability_density(n,eigen_vectors)

        # determine color based on current n value

        axs[n].plot(x_values,prob,color=(0,0.9*n/n_tot,0.8*n/n_tot))

        axs[n].set(xlabel=r'x ($\sqrt{\frac{\hbar}{m \omega}})$',ylabel=r'$|\psi(x)|^2$',title=f'n = {n}')
    

# input r values, l value, ouput: columb potential array

def coloumb_potential(r: np.ndarray,l: np.ndarray) -> np.ndarray:
    
    v = - 2/r + l*(l+1) / (r)**2

    return v


# Helpful to display this potential for a few l values

def plot_columb_potential(l_values: list, maxR: float, N: int) -> None:

    r = np.linspace(0.01,maxR,N)

    fig, ax = plt.subplots()

    ax.set(ylim=(-5,5))

    for l in l_values:

        v = coloumb_potential(r,l)

        ax.set(xlabel='$r$',ylabel='$V(r)$',title="Coloumb Potential")

        ax.plot(r,v,label=f'L = {l}')

    ax.legend()


    # input: l, max R, N points, energy eigen values and vectors

def atomic_eigen_val_vec(l: int, max_R: float, N: int) -> tuple:
    
    x_values = np.linspace(1/N,max_R,N)
    
    v_values = coloumb_potential(x_values,l) 

    # define hamiltonian
    
    H = hamiltonian(x_values,v_values)

    # compute eigen values and vectors using scipy function
    
    eigen_values, eigen_vectors = sp.linalg.eigh(H)

    return eigen_values, eigen_vectors

# input: value of quantum number l, number of expected eigen values to compute, output: first N expected eigen values

def expected_atomic_eigen_val(l: int, N_values: int) -> np.ndarray:

    # l = 0 case considered separately due to n_min = 0, hence 1/(0**2) -> inf

    if l == 0:

        expected = np.array([-np.inf])
        expected = np.append(expected,-np.array([1/(n**2) for n in range(1,N_values)]))

    else:
        expected = np.array([-1/(n**2) for n in range(l+1,N_values+(l+1))])

    return expected



# input: n_max (plots wavefunctions for n = 1, ... n_max), max radius, N points, ouput: 1D and 2D wavefunctions (2D plots are based on 1D plots - radial symmetry)

def plot_atomic_s_orbitals(n_max: int, maxR: float, N: int) -> None: 

    # create figure

    fig, axs = plt.subplots(n_max,2,figsize=(10,5*n_max))

    fig.suptitle(f'Atomic Wavefunctions, n = [1,{n_max}], l = 0 (s orbitals)')

    # define r and v

    r_values = np.linspace(1/N,maxR,N)

    v_values = coloumb_potential(r_values,0)

    # loop through each n value, starting at n = 1 (ignore n = 0 case with dirac delta at 0)

    for n in range(1,n_max+1):

        H = hamiltonian(r_values,v_values)

        eigen_val, eigen_vec = sp.linalg.eigh(H)

        # define 1D wavefunction

        prob_r = probability_density(n,eigen_vec)

        #Checked normalisation here
        
        #print(sum(prob_r))

        # define 2D wavefunction based on 1D function

        prob_x_y = np.zeros((N,N))

        # downside of 2D plot: large N can take a fair bit of time here. 

        for i in range(N):
            for j in range(N):

                # this technique isn't perfect; discussed in Part 5

                x = int(i - N/2)
                y = int(j - N/2)

                r = int(np.sqrt(x**2 + y**2))
            
                prob_x_y[i][j] = prob_r[r]
                
        # plot 1D

        lim = 10 + 20*n  # set limit based on n value

        axs[n-1][0].plot(r_values,prob_r)
        axs[n-1][0].set(xlim=(0,lim),title=f'n = {n}',ylabel='$|\psi{}(r)|^2$',xlabel='$r(a_0)$')

        # plot 2D

        axs[n-1][1].imshow(prob_x_y,extent=[-maxR,maxR,-maxR,maxR],cmap=plt.cm.YlGnBu_r)
        axs[n-1][1].set(xlim=(-lim,lim),ylim=(-lim,lim),title=f'{n}s',xlabel='$x(a_0)$',ylabel='$y(a_0)$')

# a complete orbital plotting function :) 

def plot_atomic_orbital(n: int,l: int,maxR: float,N: int) -> None: 

    # correct way to define possible m values, consqeuence of angular momentum operator.

    m_vals = np.linspace(-l,l,(2*l)+1)

    plots = int(len(m_vals))

    fig, axs = plt.subplots(plots,2,figsize=(12,6*plots))

    fig.suptitle(f'Atomic Wavefunctions, n = {n}, l = {l}, m = {min(m_vals),max(m_vals)}')

    # define r values and find probability distribution.

    r_values = np.linspace(1/N,maxR,N)
    v_values = coloumb_potential(r_values,l)


    H = hamiltonian(r_values,v_values)

    eigen_val, eigen_vec = sp.linalg.eigh(H)

    index = n - (l+1)

    # error catch if a value l+1 > n

    if index < 0:

        return 'l must be smaller than n'

    # define 1D probability density

    prob_r = probability_density(index,eigen_vec)

    #print(sum(prob_r))

    # define probability intensity in 2D x-y intensity array

    prob_x_y = np.zeros((N,N))

    for m in m_vals:

        for i in range(N):
            for j in range(N):

                x = int(i - N/2)
                y = int(j - N/2)

                # caculate theta

                theta = np.arctan2(y, x)

                # change range from [-pi,pi] to [0,2pi]:

                if theta < 0:
                    theta += 2*np.pi

                Y_l_m = sps.sph_harm(m,l,np.pi,theta) #phi = constant (can be anything) as we are only considering probality ( |phase factor|^2 = 1 )

                # compute |Y_l_m|^2

                Y_l_m_squared = np.conj(Y_l_m) * Y_l_m

                # define 2D probability based on value of r. 
                
                r = int(np.sqrt(x**2 + y**2))
                
                prob_x_y[i][j] = prob_r[r]  * np.abs(Y_l_m_squared)

        # ploting:
        
        if plots != 1:
            axs[int(m+1)][1].imshow(prob_x_y,extent=[-maxR,maxR,-maxR,maxR],cmap=plt.cm.YlGnBu_r)
            axs[int(m+1)][1].set(title=f'$(n, l, m) = ({n}, {l}, {int(m)})$',xlabel='$x(a_0)$',ylabel='$y(a_0)$')

            axs[int(m+1)][0].plot(r_values,prob_r)
            axs[int(m+1)][0].set(title=f'n = {n},l = {l}',ylabel='$|\psi{}(r)|^2$',xlabel='$r(a_0)$')
            

        else:
            axs.imshow(prob_x_y,extent=[-maxR,maxR,-maxR,maxR],cmap=plt.cm.YlGnBu_r)
            axs.set(title=f'$(n, l, m) = ({n}, {l}, {int(m)})$',xlabel='$x(a_0)$',ylabel='$y(a_0)$')


def hamiltonian_diag(x: np.ndarray,v: np.ndarray) -> tuple:

    N = len(x)

    dx = x[1] - x[0]

    on_diag = np.zeros(N)
    off_diag = np.zeros(N-1)

    for i in range(N):

        on_diag[i] = 2/((dx)**2)

        if i-1 >= 0:
            off_diag[i-1] = -1/((dx)**2)

    
    on_diag += v
        
    return on_diag, off_diag


# input: array of test N values, output: plot

def N_vs_accuracy(N_test_values: np.ndarray) -> None:

    error_list = []

    for N in N_test_values:

        x_values = np.linspace(-6,6,int(N)) 
    
        v_values = potential_qho(x_values)

        on,off = hamiltonian_diag(x_values,v_values)

        expected_vals = expected_harmonic_eigen_val(10)

        # using hamilitionian in combination with sp.linalg.eigh_tridiagonal

        eigen_values, eigen_vectors = sp.linalg.eigh_tridiagonal(on,off)

        percentage_error = abs((1 - (eigen_values[0:10]/expected_vals)) * 100)

        error_list.append(np.mean(percentage_error))

    # plot

    fig,ax = plt.subplots()

    ax.plot(N_test_values,error_list,c='red')

    ax.set(title=r'First 10 eigen values Harmonic Oscillator: N vs % Error',xlabel='N',ylabel='Error (%)')
