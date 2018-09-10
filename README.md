# delteX
delteX is a Python/Numpy script to generate an excited state structure from 
deltas, frequencies, normal modes and ground state structure.
It can also project the dimensionless displacements onto internal coordinates with dimensions.

To install:

pip install delteX/

Then to run, need the following files in the working directory:

normal_modes.dat # Normal mode displacement vectors
freqs.dat # normal mode frequencies
shifts_dnc.dat # Dimensionless normal mode displacements
xyz.dat # cartesian coordinate structure of molecule
masses.dat # atomic mass of each atom 

See example data for the required structure of each file. Eventually file parsing scripts for ORCA output will be incorporated directly into delteX. This is on the to do list.

With all files in the same directory and delteX installed, run:

delteX a b c d

Where "a", "b", "c", and "d" are the coordinate numbers of the atoms which define the internal coordinate of interest.

The program will output the total displacement along all possible internal coordinates (dihedrals, angles, and bonds) that can be created with this set of 4 atoms. Additionally, separate files will be created which contain the internal coordinate displacement projected onto normal coordinates which can be used to visualize the contribution to displcament along a specific coordinate from different peaks in a Raman spectrum.

The program will also output the Franck-Condon relaxed excited state structure determined by projecting the dimensionless normal mode displacements onto cartesian coordinates. This is the excited state geometry as predicted under the Independent Mode Displaced Harmonic Oscillator (IMDHO) model.
