import numpy as np
from math import *
from numpy.linalg import inv

class wilson_B:

	def __init(self):
		pass

## Load in the required data extracted from frequency and es gradient calculation

	def mass(self,N):	
	mass = np.asmatrix(np.loadtxt('masses.dat'))
		mass_matrix = np.zeros((3*N,3*N))
		for i in range(N):
			j = 3*i
			mass_matrix[j][j] = (np.sqrt(np.transpose(mass)[i]))
			mass_matrix[j+1][j+1] = (np.sqrt(np.transpose(mass)[i]))
			mass_matrix[j+2][j+2] = (np.sqrt(np.transpose(mass)[i]))
		return mass_matrix

	def norm_modes(self,N):
		n_m =  self.mass(N)*np.asmatrix(np.loadtxt('normal_modes.dat',usecols=range(0,51)))
		return np.asmatrix(n_m)

	def shifts_dnc(self):
		shifts = np.zeros(51)
		shifts_dnc = np.loadtxt('shifts_dnc.dat')
		for i in range(51):
			shifts[i] = shifts_dnc[i]
		return shifts

	def freqs(self):
		freq = np.zeros(51)
		freqs = np.array(np.loadtxt("simFreq.dat")) 
		for i in range(51):
			freq[i] = freqs[i] 
		for i in range(len(freq)):
			if freq[i]!=0.0:
				freq[i] = np.sqrt(1./freq[i])
		return freq

	def xyz(self,N):
		xyz_raw = np.asmatrix(np.loadtxt('xyz.dat',usecols=range(3)))*0.52917724900001 # 0.52  converts bohr to angstroms
		xyz = np.zeros((3*N,1))
		for i in range(N):
			xyz[3*i] = np.transpose(xyz_raw[i])[0]
			xyz[3*i+1] = np.transpose(xyz_raw[i])[1]
			xyz[3*i+2] = np.transpose(xyz_raw[i])[2]
		return xyz

	def reform_xyz(self,xyz_col,N):
		xyz = np.zeros((N,3))
		xyz_col = np.asarray(xyz_col)
		for j in range(N):
			xyz[j][0] = xyz_col[0][3*j]
			xyz[j][1] = xyz_col[0][3*j+1]
			xyz[j][2] = xyz_col[0][3*j+2]
		return xyz
	
## Calculate the dihedral angle of the input 4 atoms given in Cartesian coordinates ##

	# Unit vectors
	def e12(self,coords):
		return (coords[1,:] - coords[0,:])/np.sqrt(np.abs(coords[1,0] - coords[0,0])**2+np.abs(coords[1,1] - coords[0,1])**2+np.abs(coords[1,2] - coords[0,2])**2)
	def e23(self,coords):
		return (coords[2,:] - coords[1,:])/np.sqrt(np.abs(coords[2,0] - coords[1,0])**2+np.abs(coords[2,1] - coords[1,1])**2+np.abs(coords[2,2] - coords[1,2])**2)
	def e34(self,coords):
		return (coords[3,:] - coords[2,:])/np.sqrt(np.abs(coords[3,0] - coords[2,0])**2+np.abs(coords[3,1] - coords[2,1])**2+np.abs(coords[3,2] - coords[2,2])**2)
	# Angles 
	def phi2(self,e12,e23):
		return np.arccos(np.dot(e12,e23))
	def phi3(self,e23,e34):
		return np.arccos(np.dot(e23,e34))
	# Calculate Dihedral
	def B_dihedral(self,coords):
		return (180.0/np.pi)*np.arccos(np.dot(np.cross(self.e12(coords),self.e23(coords)),np.cross(self.e23(coords),self.e34(coords)))/(np.sin(self.phi2(self.e12(coords),self.e23(coords)))*np.sin(self.phi3(self.e23(coords),self.e34(coords)))))

## Calculate Bond Distances given two atoms in Cartesian Coordinates ##

	def bond23(self,coords):
		return (coords[2,:] - coords[1,:])

####################################################################################

## Dihedralize ##
## Takes 4 atom numbers as input and the number of modes total ##
###################################################################################
	def internal_matrix(self,matrix,ax,bx,cx,dx,N,internal):

		## Next dihedralize the inverted Hessian and energy gradient ##
		matrix_dihX = np.zeros((3*N,4,3)) #3*N
		matrix_int = np.zeros(3*N) #3*N
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
					matrix_dihX[j,0,0] = matrix[i,j] 
				if i == by:
					matrix_dihX[j,1,1] = matrix[i,j]
				if i == cz:
					matrix_dihX[j,2,2] = matrix[i,j]
				if i == dx:
					matrix_dihX[j,3,0] = matrix[i,j]
				if i == ay:
					matrix_dihX[j,0,1] = matrix[i,j]
				if i == bz:
					matrix_dihX[j,1,2] = matrix[i,j] 
				if i == cx:
					matrix_dihX[j,2,0] = matrix[i,j] 
				if i == dy:
					matrix_dihX[j,3,1] = matrix[i,j] 
				if i == az:
					matrix_dihX[j,0,2] = matrix[i,j]
				if i == bx:
					matrix_dihX[j,1,0] = matrix[i,j] 
				if i == cy:
					matrix_dihX[j,2,1] = matrix[i,j]
				if i == dz:
					matrix_dihX[j,3,2] = matrix[i,j]
			if internal == "dihedral":
				matrix_int[j] = self.B_dihedral(matrix_dihX[j,:,:])
			elif internal == "bond23":
				matrix_int[j] = np.sqrt(self.bond23(matrix_dihX[j,:,:])[0]**2+self.bond23(matrix_dihX[j,:,:])[1]**2+self.bond23(matrix_dihX[j,:,:])[2]**2)
		return matrix_int

	def internal_array(self,array,ax,bx,cx,dx,N,internal):
		array_dihX = np.zeros((4,3))
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
				array_dihX[0,0] = array[i,0]
			if i == by:
				array_dihX[1,1] = array[i,0]
			if i == cz:
				array_dihX[2,2] = array[i,0]
			if i == dx:
				array_dihX[3,0] = array[i,0]
			if i == ay:
				array_dihX[0,1] = array[i,0]
			if i == bz:
				array_dihX[1,2] = array[i,0]
			if i == cx:
				array_dihX[2,0] = array[i,0]
			if i == dy:
				array_dihX[3,1] = array[i,0]
			if i == az:
				array_dihX[0,2] = array[i,0]
			if i == bx:
				array_dihX[1,0] = array[i,0]
			if i == cy:
				array_dihX[2,1] = array[i,0]
			if i == dz:
				array_dihX[3,2] = array[i,0]
		if internal == "dihedral":
			array_int = self.B_dihedral(array_dihX)
		elif internal == "bond23":
			array_int = np.sqrt(self.bond23(array_dihX)[0]**2+self.bond23(array_dihX)[1]**2+self.bond23(array_dihX)[2]**2)
		return array_int
	
	def delta(self,eigen_hess,norm_modes,mass,en_grad,N):
		delta_int = np.zeros(3*N)
		#print np.shape(eigen_hess), np.shape(norm_modes), np.shape(mass), np.shape(en_grad)
		#exit()
		delta_int = eigen_hess*np.transpose(norm_modes)*(1/(mass**0.5))*en_grad
		return delta_int
	
	def deltaX2deltaInt(self,xyz,mass,shifts_mwnc,ax,bx,cx,dx,N):
		shifts_nc = np.zeros(N)
		shifts_nc = mass*(xyz+shifts_mwnc)
		return self.dihedral_array(shifts_nc,ax,bx,cx,dx,N)


	# Convert cartesian/normal displacement matrix to internal/normal
	def n_m_CN2IN_dihed(self,n_m,ax,bx,cx,dx,N):	
		return self.internal_matrix(n_m,ax,bx,cx,dx,N,"dihedral")

	def n_m_CN2IN_bond(self,n_m,ax,bx,cx,dx,N):
		return self.internal_matrix(n_m,ax,bx,cx,dx,N,"bond23")

	# Gives the dihedral contribution 	
	def deltaExp2deltaInt_dihed(self,shifts_dnc,n_m,ax,bx,cx,dx,N,xyz):
		np.savetxt("A_tao.dat",self.n_m_CN2IN_dihed((n_m+xyz),ax,bx,cx,dx,N))
		xyz_equil=self.internal_array(xyz,ax,bx,cx,dx,N,"dihedral") 
		return 5.8065*np.dot((self.n_m_CN2IN_dihed((n_m+xyz),ax,bx,cx,dx,N)-xyz_equil),self.freqs()*shifts_dnc)
	
	def deltaExp2deltaInt_bond(self, shifts_dnc,n_m,ax,bx,cx,dx,N,xyz):
		np.savetxt("A_taoBond.dat",self.n_m_CN2IN_bond((n_m),ax,bx,cx,dx,N))
		xyz_equil = self.internal_array(xyz,ax,bx,cx,dx,N,"bond23")
		return 5.8065*np.dot((self.n_m_CN2IN_bond((n_m+xyz),ax,bx,cx,dx,N)-xyz_equil),self.freqs()*shifts_dnc)

	def delta2xyzESgeom(self,n_m,shifts_dnc):
		return 5.8065*np.dot(n_m,self.freqs()*shifts_dnc)


	def force2forceInt(self,force,ax,bx,cx,dx,N):
		return self.dihedral_array(force,ax,bx,cx,dx,N)

	
	def action_list(self,ax,bx,cx,dx,N):
		## Initialize the required data ##
		norm_modes = self.norm_modes(N)
		xyz = self.xyz(N)
		shifts_dnc = self.shifts_dnc()	

		## Compute equilibrium value for internal coordinate ##
		equil_dihed = self.internal_array(xyz,ax,bx,cx,dx,N,"dihedral")
		equil_bond23 = self.internal_array(xyz,ax,bx,cx,dx,N,"bond23")

		## Compute excited state shift in xyz ##
		es_struc = self.delta2xyzESgeom(norm_modes,shifts_dnc)
		
		## Reform xyz es structure for visualisation ##
		es_xyz_reform = self.reform_xyz((es_struc+np.transpose(xyz)),N)
		gs_xyz_reform = self.reform_xyz((np.transpose(xyz)),N)
		np.savetxt("es_struc.dat",es_xyz_reform)
		np.savetxt("gs_struc.dat",gs_xyz_reform)	
	
		print self.deltaExp2deltaInt_dihed(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz) 
		print self.deltaExp2deltaInt_bond(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz) 

		print equil_dihed
		print equil_bond23


wilson_B = wilson_B()
print "DMA"	
wilson_B.action_list(16,15,10,11,46) 
print "Ar-Core"
wilson_B.action_list(1,0,38,39,46)

# O1 Ar-Core  1,0,22,26 #DMA 12,13,14,18 49atoms
# O2 Ar-Core 1,0,38,39 # DMA 16,15,10,11 46atoms





