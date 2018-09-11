import numpy as np
from math import *
from numpy.linalg import inv
import warnings

warnings.filterwarnings('ignore') # Supresses 'casting to real discards complex part' warning

# masses of elements #
def mass(N):	
	mass = np.asmatrix(np.loadtxt('masses.dat'))
	mass_matrix = np.zeros((3*N,3*N))
	for i in range(N):
		j = 3*i
		mass_matrix[j][j] = (np.sqrt(np.transpose(mass)[i]))
		mass_matrix[j+1][j+1] = (np.sqrt(np.transpose(mass)[i]))
		mass_matrix[j+2][j+2] = (np.sqrt(np.transpose(mass)[i]))
	return mass_matrix

# normal mode displacement matrix #
def norm_modes(N):
	n_m =  mass(N)*np.asmatrix(np.loadtxt('normal_modes.dat')) 
	return np.asmatrix(n_m)

# dimensionless normal mode displacements
def shifts_dnc(N):
	shifts = np.zeros(3*N) # Replace 3*N with specific number if using a custom set of normal modes
	shifts_dnc = np.loadtxt('shifts_dnc.dat')
	for i in range(3*N): # Replace 3*N with specific number if using a custom set of normal modes
		shifts[i] = shifts_dnc[i]
	return shifts

def freqs(N):
	freq = np.zeros(3*N) # Replace 3*N with specific number if using a custom set of normal modes
	freqs = np.array(np.loadtxt("freqs.dat")) 
	for i in range(3*N): # Replace 3*N with specific number if using a custom set of normal modes
		freq[i] = freqs[i] 
	for i in range(len(freq)):
		if freq[i]!=0.0:
			freq[i] = np.sqrt(1./freq[i])
	return freq

def xyz(N):
	xyz_raw = np.asmatrix(np.loadtxt('xyz.dat',usecols=range(3)))*0.52917724900001 # 0.52  converts bohr to angstroms
	xyz = np.zeros((3*N,1))
	for i in range(N):
		xyz[3*i] = np.transpose(xyz_raw[i])[0]
		xyz[3*i+1] = np.transpose(xyz_raw[i])[1]
		xyz[3*i+2] = np.transpose(xyz_raw[i])[2]
	return xyz

def reform_xyz(xyz_col,N):
	xyz = np.zeros((N,3))
	xyz_col = np.asarray(xyz_col)
	for j in range(N):
		xyz[j][0] = xyz_col[0][3*j]
		xyz[j][1] = xyz_col[0][3*j+1]
		xyz[j][2] = xyz_col[0][3*j+2]
	return xyz

## Calculate the internal coordinates of the input 4 atoms given in Cartesian coordinates ##

# Unit vectors
def e12(coords):
	return (coords[1,:] - coords[0,:])/np.sqrt(np.abs(coords[1,0] - coords[0,0])**2+np.abs(coords[1,1] - coords[0,1])**2+np.abs(coords[1,2] - coords[0,2])**2)
def e23(coords):
	return (coords[2,:] - coords[1,:])/np.sqrt(np.abs(coords[2,0] - coords[1,0])**2+np.abs(coords[2,1] - coords[1,1])**2+np.abs(coords[2,2] - coords[1,2])**2)
def e34(coords):
	return (coords[3,:] - coords[2,:])/np.sqrt(np.abs(coords[3,0] - coords[2,0])**2+np.abs(coords[3,1] - coords[2,1])**2+np.abs(coords[3,2] - coords[2,2])**2)
# Angle 
def phi(e1,e2):
	return np.arccos(np.dot(e1,e2))
# Dihedral
def B_dihedral(coords):
	return (180.0/np.pi)*np.arccos(np.dot(np.cross(e12(coords),e23(coords)),np.cross(e23(coords),e34(coords)))/(np.sin(phi(e12(coords),e23(coords)))*np.sin(phi(e23(coords),e34(coords)))))
# Bond distances

def bond23(coords):
	return (coords[2,:] - coords[1,:])
def bond12(coords):
	return (coords[1,:] - coords[0,:])
def bond34(coords):
	return (coords[3,:] - coords[2,:])

# Form matrix elements of coordinates #
def internal_matrix(matrix,ax,bx,cx,dx,N,internal):

	matrix_X = np.zeros((3*N,4,3)) 
	matrix_int = np.zeros(3*N) 
	ax=ax*3
	ay=ax+1
	az=ax+2
	bx=bx*3
	by=bx+1
	bz=bx+2
	cx=cx*3
	cy=cx+1
	cz=cx+2
	dx=dx*3
	dy=dx+1
	dz=dx+2
	for j in range(24):
		for i in range(np.size(matrix,0)):
			if i == ax:
				matrix_X[j,0,0] = matrix[i,j] 
			if i == by:
				matrix_X[j,1,1] = matrix[i,j]
			if i == cz:
				matrix_X[j,2,2] = matrix[i,j]
			if i == dx:
				matrix_X[j,3,0] = matrix[i,j]
			if i == ay:
				matrix_X[j,0,1] = matrix[i,j]
			if i == bz:
				matrix_X[j,1,2] = matrix[i,j] 
			if i == cx:
				matrix_X[j,2,0] = matrix[i,j] 
			if i == dy:
				matrix_X[j,3,1] = matrix[i,j] 
			if i == az:
				matrix_X[j,0,2] = matrix[i,j]
			if i == bx:
				matrix_X[j,1,0] = matrix[i,j] 
			if i == cy:
				matrix_X[j,2,1] = matrix[i,j]
			if i == dz:
				matrix_X[j,3,2] = matrix[i,j]
		if internal == "dihedral":
			matrix_int[j] = B_dihedral(matrix_X[j,:,:])
		elif internal == "angle123":
			matrix_int[j] = phi(e12(matrix_X[j,:,:]),e23(matrix_X[j,:,:]))
		elif internal == "angle234":
			matrix_int[j] = phi(e23(matrix_X[j,:,:]),e34(matrix_X[j,:,:]))
		elif internal == "bond23":
			matrix_int[j] = np.sqrt(bond23(matrix_X[j,:,:])[0]**2+bond23(matrix_X[j,:,:])[1]**2+bond23(matrix_X[j,:,:])[2]**2)
		elif internal == "bond12":
			matrix_int[j] = np.sqrt(bond12(matrix_X[j,:,:])[0]**2+bond12(matrix_X[j,:,:])[1]**2+bond12(matrix_X[j,:,:])[2]**2)
		elif internal == "bond34":
			matrix_int[j] = np.sqrt(bond34(matrix_X[j,:,:])[0]**2+bond34(matrix_X[j,:,:])[1]**2+bond34(matrix_X[j,:,:])[2]**2)
	return matrix_int

# form array elements of coordinates #
def internal_array(array,ax,bx,cx,dx,N,internal):
	array_X = np.zeros((4,3))
	ax=ax*3
	ay=ax+1
	az=ax+2
	bx=bx*3
	by=bx+1
	bz=bx+2
	cx=cx*3
	cy=cx+1
	cz=cx+2
	dx=dx*3
	dy=dx+1
	dz=dx+2
	for i in range(np.size(array,0)):
		if i == ax:
			array_X[0,0] = array[i,0]
		if i == by:
			array_X[1,1] = array[i,0]
		if i == cz:
			array_X[2,2] = array[i,0]
		if i == dx:
			array_X[3,0] = array[i,0]
		if i == ay:
			array_X[0,1] = array[i,0]
		if i == bz:
			array_X[1,2] = array[i,0]
		if i == cx:
			array_X[2,0] = array[i,0]
		if i == dy:
			array_X[3,1] = array[i,0]
		if i == az:
			array_X[0,2] = array[i,0]
		if i == bx:
			array_X[1,0] = array[i,0]
		if i == cy:
			array_X[2,1] = array[i,0]
		if i == dz:
			array_X[3,2] = array[i,0]
	if internal == "dihedral":
		array_int = B_dihedral(array_X)
	elif internal == "angle123":
		array_int = phi(e12(array_X),e23(array_X))
	elif internal == "angle234":
		array_int = phi(e23(array_X),e34(array_X))
	elif internal == "bond12":
		array_int = np.sqrt(bond12(array_X)[0]**2+bond12(array_X)[1]**2+bond12(array_X)[2]**2)
	elif internal == "bond23":
		array_int = np.sqrt(bond23(array_X)[0]**2+bond23(array_X)[1]**2+bond23(array_X)[2]**2)
	elif internal == "bond34":
		array_int = np.sqrt(bond34(array_X)[0]**2+bond34(array_X)[1]**2+bond34(array_X)[2]**2)
	return array_int


# Convert cartesian/normal displacement matrix to internal/normal
def n_m_CN2IN_dihed(n_m,ax,bx,cx,dx,N):	
	return internal_matrix(n_m,ax,bx,cx,dx,N,"dihedral") #matrix

def n_m_CN2IN_angle(n_m,ax,bx,cx,dx,N,angle):	
	return internal_matrix(n_m,ax,bx,cx,dx,N,angle)

def n_m_CN2IN_bond(n_m,ax,bx,cx,dx,N,bond):
	return internal_matrix(n_m,ax,bx,cx,dx,N,bond)

def delta2xyzESgeom(n_m,shifts_dnc,N):
	return 5.8065*np.dot(n_m,freqs(N)*shifts_dnc)

def dihed_shifts_2(shifts_dnc,n_m,ax,bx,cx,dx,N,xyz,equil_dihed):
	trans = np.transpose(np.matrix(5.8065*np.multiply(n_m,np.transpose(freqs(N)*shifts_dnc))))
	return n_m_CN2IN_dihed(np.transpose(trans)+xyz,ax,bx,cx,dx,N)-equil_dihed

def bond_shifts_2(shifts_dnc,n_m,ax,bx,cx,dx,N,xyz,equil_bond,bond):
	trans = np.matrix(5.8065*np.multiply(n_m,np.transpose(freqs(N)*shifts_dnc)))
	return n_m_CN2IN_bond(trans+xyz,ax,bx,cx,dx,N,bond)#-equil_bond


def dihed_shifts(shifts_dnc,n_m,ax,bx,cx,dx,N,xyz,equil_dihed):
	return 5.8065*np.transpose(n_m_CN2IN_dihed((n_m+xyz),ax,bx,cx,dx,N))*freqs(N)*shifts_dnc

def angle_shifts(shifts_dnc,n_m,ax,bx,cx,dx,N,xyz,equil_angle,angle):
	return 5.8065*np.transpose(n_m_CN2IN_angle((n_m+xyz),ax,bx,cx,dx,N,angle))*freqs(N)*shifts_dnc

def bond_shifts(shifts_dnc,n_m,ax,bx,cx,dx,N,xyz,equil_bond,bond):
	return 5.8065*np.transpose(n_m_CN2IN_bond((n_m+xyz),ax,bx,cx,dx,N,bond)-equil_bond)*freqs(N)*shifts_dnc


	



