#biblioteki
from fenics import *
import numpy as np

def coor(line):
	line_s = line.split('\t')
	return [float(line_s[0]), float(line_s[1]), float(line_s[2])] 

def tab_vessel(name):

	#result
	result = []

	#mesh settings - possible to change
	nx = 20
	ny = 20
	nz = 20

	mesh = 	BoxMesh(Point(0,0,0), Point(6,6,6), nx, ny, nz)
	
	P1 = FiniteElement('P', tetrahedron, 1)
	V = FunctionSpace(mesh, P1)

	#function we want to define
	k = Function(V)

	#load
	with open(name+".txt") as f:
		data =  (f.read().splitlines())
	f.close()

	#for each cell
	for cell_no in range(mesh.num_cells()):

		print 100*float(cell_no)/mesh.num_cells() 

		assigned = 0

		#for each element check if cell contains point
		for i in range(0, len(data)):

			coordinates = coor(data[i])
			pnt = Point(coordinates[0], coordinates[1], coordinates[2])

			if Cell(mesh, cell_no).contains(pnt):
				result.append(1)
				assigned = 1
				break

		#if not assigned
		if assigned == 0:
			result.append(0)
	

	#save
	w = open("f_expression_ver" + name + ".txt", "a")
	for i in range(0, len(result)):
		w.write(str(result[i]))

	w.close()

#here put name of the file	
tab_vessel("sector_kol_7_13_3_9_3.5_9.5")