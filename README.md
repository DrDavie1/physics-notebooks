# Quantum Eigenstates ψ

This problem was set as a university assignment within a jupyter notebook, however, I have removed specific assignement details from the notebook and will give details here. The main_code.py file contains the functions used to solve the problem. 

References:

- [Griffiths, D.J. and Schroeter, D.F., 2018. Introduction to quantum mechanics. Cambridge university press](https://www.google.co.uk/books/edition/Introduction_to_Quantum_Mechanics/LWRnDwAAQBAJ?hl=en&gbpv=0).
  - Chapter 2.3: The Harmonic Oscillator
  - Chapter 4.2: The Hydrogen Atom

#

The code consisted of 4 parts including an extension.

### Part 1

Plots a harmonic pontenial defined by

$$ V(x) = x^2 $$

### Part 2

Defines a Hamiltonian matrix function. 

The hamiltonian is given by: 

$$ H = \frac{-\hbar^2}{2m}\frac{d^2}{dx^2} + V $$

Applying the finite difference method and using dimensionless units, the KE matrix of the hamilonian $D$ is defined by:

$$D_{i,i} = \frac{2}{(\Delta x)^2}$$

$$D_{i,i+1} = D_{i,i-1} = \frac{-1}{(\Delta x)^2}$$

Hence:

$$ H = D+V $$


### Part 3

Finds the eigen values and states for a Quantum Harmonic oscillator in 1D.

<img src="https://github.com/DrDavie1/quantum-eigenstates/blob/main/Media/harmosc.png" width="50%" height="50%">


### Part 4

Finds the eigen values and states of a hydrogen atom with l = 0 in 1D and plots in 2D.

#

<img src="https://github.com/DrDavie1/quantum-eigenstates/blob/main/Media/3s.png" width="40%" height="40%">

### Extension

- Finds the eigen values and states of a hydrogen atom with l > 0

#

<img src="https://github.com/DrDavie1/quantum-eigenstates/blob/main/Media/p.png" width="40%" height="40%">

- Compares the efficiency of scipy eigen problem solvers.
    
