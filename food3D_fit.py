#biblioteki
from pyswarm import pso
from fenics import * 
import numpy as np
import time

def rd_solver(A, e, M, ptsx, ptsy, ptsz, typ):

	#cost function
	cost = 0	

	#data - time
	T = 90.0			#total time of simulation
	num_steps = 180		#number of steps
	dt = T/num_steps	#one time step 

	#data
	d_n = Constant(0.56)		#dyfuzja dla skladnikow odzywcz
	k = Constant(dt)		#sie przyda

	#mesh, element...
	print "Loading function space, geometry ..."
	nx = 20
	ny = 20
	nz = 20
	mesh = 	BoxMesh(Point(0,0,0), Point(6,6,6), nx, ny, nz)
	P1 = FiniteElement('P', tetrahedron, 1)
	V = FunctionSpace(mesh, P1)
	
	#test functions and solution
	v = TestFunction(V)
	u = Function(V)
	
	#initial conditions
	print "Loading initial conditions ..."
	u_n = Function(V)
	u_0 = Expression('1.17', degree=1)
	u_n = project(u_0, V)

	#load vessels data
	f = CellFunction("int", mesh, 3)

	if typ == "r":
		with open("f_r.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()

	elif typ == "g":
		with open("f_g.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
	
	elif typ == "b":
		with open("f_b.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
	r=0

	for cell in cells(mesh):
		f[r] = int(data[r])
		r = r+1

	class MyExpression1(Expression):

		def __init__(self, f, **kwargs):
			self.f = f
		
		def eval_cell(self, value, x, ufc_cell):
			value[0] = self.f[ufc_cell.index]

	f0 = MyExpression1(f, degree=5)

	#variational formulation 
	F = ((u - u_n)*v/k)*dx + d_n*dot(grad(u), grad(v))*dx + e*u*v*(((u*u))/(M+((u*u))))*dx-A*f0*v*dx

	#picture saving (disabled)
	#vtk = File('cancer3D/food_density.pvd')

	#progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)
	
	#solve
	t = 0
	for n in range(num_steps):
		
		# update time
		t += dt
	
		# solve
		solve(F == 0, u)
	
		# save picture (disabled)
		#vtk << (u, t)

		# compute cost in supremum norm
	
		for i in range(0, len(ptsx)):
			
			cost = max((u(ptsx[i], ptsy[i], ptsz[i]) - 1.17)**2, cost)

		# aktualizacja u_n, u_n1, u_n2 i paska
		u_n.assign(u)
		progress.update(t / T)

	return cost


def cost_function(x):

	#parameters to fit
	e = x[0]
	M = x[1]
	Ar = x[2]
	Ag = x[3]
	Ab = x[4]

	#at this points we evaluate cost function
	#for red
	r_x=[3,1,3,5,3,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,2.5,2.5,2.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,2.5,2.5,2.5]
	r_y=[3,3,3,3,3,1.75,5,5,5,1,1,1,2,2,2,2,2,2,4,4,4,5,5,5,1,3,5,1,3,5,1,3,5,4,4,4,4,4,4]
	r_z=[3,3,5,5,2,3,1,1,1,5,5,5,5,5,5,4,4,4,1,1,1,2,2,2,0.8,0.8,0.8,2.5,2.5,2.5,5,5,5,1,3,5,1,3,5]
	
	#for green
	g_x=[5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
	g_y=[1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,3,4,5,2,3,4,5,1.3,2,3,4,5,1,2,3,4,5,1,1,1,1,3,3,3,3,5,5,5,5,1,1]
	g_z=[1,1,1,3,3,3,5,5,5,1,1,1,3,3,3,5,5,5,2,2,2,3,3,3,5,5,5,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5,2,3,4,5,2,3,4,5,2,3,4,5,1,1]

	#for blue
	b_x=[2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,4]
	b_y=[1,1,1,1,3,3,3,3,5,5,5,5,1,2,3,4,5,1,2,3,4,5,1,1,1,1,3,3,3,3,5,5,5,5,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,1]
	b_z=[1,2,3,4,1,2,3,4,1,2,3,4,4,4,4,4,4,1,1,1,1,1,1.25,3,3.5,4.5,1.25,3,3.5,4.5,1.25,3,3.5,4.5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1.5]

	cost = 0
	cost = cost + rd_solver(Ar, e, M, r_x,r_y, r_z,"r") + rd_solver(Ag, e, M, g_x, g_y, g_z,"g") + rd_solver(Ab, e, M, b_x, b_y, b_z, "b")
	w = open("results.txt","a")
	w.write(str(e) + "\t" + str(M) + "\t" + str(Ar) + "\t" + str(Ag) + "\t" + str(Ab) +  "\t" + str(cost) + "\n")
	w.close()
	return cost

xopt, fopt = pso(cost_function, [0.04,0.82,1.75,2.40,2.0], [0.10,0.95,1.90,2.56,2.16], debug = True)