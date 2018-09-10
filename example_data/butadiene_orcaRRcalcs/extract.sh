#!/bin/bash

FILE_HESS="buta_tddft.rrhess"
NATOMS=10

# Extract Cartesian Displacements

for i in $(seq 0 5); do
	grep -F '^n^' -A 756 $FILE_HESS | sed -n -e $[i*31+2],$[i*31+31]p | awk '{print $2,$3,$4,$5,$6,$7}' > modes_$i.tmp
done 

paste modes_{0..5}.tmp | column -s $'\t' -t  > normal_modes.dat

rm *modes_*

# Extract Vibrational frequencies
grep '$vibrational_frequencies' -A 33 $FILE_HESS | sed -n -e 3,33p | awk '{print $2}' > freqs.dat 

# Extract Hessian

#for i in $(seq 0 24); do
#	grep -F '^h^' -A 3700 $FILE_HESS | sed -n -e $[i*148+2],$[i*148+148]p | awk '{print $2,$3,$4,$5,$6,$7}' > modes_$i.tmp
# done

#paste modes_{0..24}.tmp | column -s $'\t' -t  > hessian.dat

#rm *modes_*

# Extract the cartesian gradient

#grep '$energy_cartesian_gradients' -A 150 $FILE_HESS | sed -n -e 2,150p | awk '{print $3}' > es_gradient.dat

# Extract Eigenvalues of the mass weighted hessian

#grep '$eigenvalues_mass_weighted_hessian' -A 150 $FILE_HESS | sed -n -e 3,150p | awk '{print $2}' > eigen_hess.dat 

# Extract masses 

grep '$atoms' -A 13 $FILE_HESS | sed -n -e 3,34p | awk '{print $2}' > masses.dat 

# Extract elements

grep '$atoms' -A 13 $FILE_HESS | sed -n -e 3,34p | awk '{print $1}' > elements.dat 
#paste elements.dat es_struc.dat > O1_es_struc.dat
#paste elements.dat gs_struc.dat > O1_gs_struc.dat

# Extract xyz 

grep '$atoms' -A 33 $FILE_HESS | sed -n -e 3,34p | awk '{print $3,$4,$5}' > xyz.dat 

# extract shifts_mwnc
grep '$shifts_dnc' -A 33 $FILE_HESS | sed -n -e 4,33p | awk '{print $2}' > shifts_dnc.dat 


