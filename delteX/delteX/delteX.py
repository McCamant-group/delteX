def main():

	import numpy as np
	from math import *
	from numpy.linalg import inv
	import sys
	import wilsonB as w

	ax = int(sys.argv[1])
	bx = int(sys.argv[2])
	cx = int(sys.argv[3])
	dx = int(sys.argv[4])
	N = int(sys.argv[5])

	## Initialize the required data ##
	norm_modes = w.norm_modes(N)
	xyz = w.xyz(N)
	shifts_dnc = w.shifts_dnc(N)	

	## Compute equilibrium value for internal coordinate ##
	equil_dihed1234 = w.internal_array(xyz,ax,bx,cx,dx,N,"dihedral1234")
	equil_angle123 = (180.0/pi)*w.internal_array(xyz,ax,bx,cx,dx,N,"angle123")
	equil_angle234 = (180.0/pi)*w.internal_array(xyz,ax,bx,cx,dx,N,"angle234")
	equil_bond12 = w.internal_array(xyz,ax,bx,cx,dx,N,"bond12")
	equil_bond23 = w.internal_array(xyz,ax,bx,cx,dx,N,"bond23")
	equil_bond34 = w.internal_array(xyz,ax,bx,cx,dx,N,"bond34")

	equil_dihedMat = w.n_m_CN2IN(norm_modes+xyz,ax,bx,cx,dx,N,"dihedral1234")-equil_dihed1234
	for i in range(len(equil_dihedMat)):
		if equil_dihedMat[i] == -1.0*equil_dihed1234:
			equil_dihedMat[i] += equil_dihed1234

	## Compute excited state shift in xyz ##
	es_struc = w.delta2xyzESgeom(norm_modes,shifts_dnc,N)
	
	## Reform xyz es structure for visualisation ##
	es_xyz = np.transpose((es_struc+np.transpose(xyz)))
	es_xyz_reform = w.reform_xyz(np.transpose(es_xyz),N)
	gs_xyz_reform = w.reform_xyz((np.transpose(xyz)),N)
	np.savetxt("es_struc.dat",es_xyz_reform)
	np.savetxt("gs_struc.dat",gs_xyz_reform)	

	es_dihed1234 = w.internal_array(es_xyz,ax,bx,cx,dx,N,"dihedral1234")
	es_angle123 = (180.0/pi)*w.internal_array(es_xyz,ax,bx,cx,dx,N,"angle123")
	es_angle234 = (180.0/pi)*w.internal_array(es_xyz,ax,bx,cx,dx,N,"angle234")
	es_bond12 = w.internal_array(es_xyz,ax,bx,cx,dx,N,"bond12")
	es_bond23 = w.internal_array(es_xyz,ax,bx,cx,dx,N,"bond23")
	es_bond34 = w.internal_array(es_xyz,ax,bx,cx,dx,N,"bond34")

	dihedral1234 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_dihed1234,"dihedral1234")		
	angle123 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_angle123,"angle123")
	angle234 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_angle234,"angle234")
	bond12 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_bond12,"bond12")
	bond23 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_bond23,"bond23")
	bond34 = w.int_shifts(shifts_dnc,norm_modes,ax,bx,cx,dx,N,xyz,equil_bond34,"bond34")


	coord1234 = str(ax)+"_"+str(bx)+"_"+str(cx)+"_"+str(dx)
	coord123 = str(ax)+"_"+str(bx)+"_"+str(cx)
	coord234 = str(bx)+"_"+str(cx)+"_"+str(dx)
	coord12 = str(ax)+"_"+str(bx)
	coord23 = str(bx)+"_"+str(cx)
	coord34 = str(cx)+"_"+str(dx)
	np.savetxt("dihedral_"+coord1234+"_GSnm.dat",equil_dihedMat)
	np.savetxt("dihedral_"+coord1234+"_nm.dat",dihedral1234)
	np.savetxt("angle_"+coord123+"_nm.dat",angle123)
	np.savetxt("angle_"+coord234+"_nm.dat",angle234)
	np.savetxt("bond_"+coord12+"_nm.dat",bond12)
	np.savetxt("bond_"+coord23+"_nm.dat",bond23)
	np.savetxt("bond_"+coord34+"_nm.dat",bond34)


	with open("deltex_"+coord1234+".out",'w') as f:
		f.write("Ground state equilibrium geometry \n")
		f.write("Dihedral "+coord1234+" = "), f.write(str(equil_dihed1234)), f.write(" degrees \n")	
		f.write("Angle "+coord123+" = "), f.write(str(equil_angle123)), f.write(" degrees \n")
		f.write("Angle "+coord234+" = "), f.write(str(equil_angle234)), f.write(" degrees \n")
		f.write("Bond "+coord12+" = "), f.write(str(equil_bond12)), f.write(" Angstroms \n")
		f.write("Bond "+coord23+" = "), f.write(str(equil_bond23)), f.write(" Angstroms \n")
		f.write("Bond" +coord34+" = "), f.write(str(equil_bond34)), f.write(" Angstroms \n\n")
		f.write("Excited State geometry \n")
		f.write("ES Dihedral "+coord1234+" = "), f.write(str(es_dihed1234)), f.write(" degrees \n")	
		f.write("ES Angle "+coord123+" = "), f.write(str(es_angle123)), f.write(" degrees \n")
		f.write("ES Angle "+coord234+" = "), f.write(str(es_angle234)), f.write(" degrees \n")
		f.write("ES Bond "+coord12+" = "), f.write(str(es_bond12)), f.write(" Angstroms \n")
		f.write("ES Bond "+coord23+" = "), f.write(str(es_bond23)), f.write(" Angstroms \n")
		f.write("ES Bond "+coord34+" = "), f.write(str(es_bond34)), f.write(" Angstroms \n")
	f.close()
