#biblioteki
import numpy as np
from fenics import * 
import time
from pyswarm import pso
import random

def rd_solver(a, g, d, d_c, typ, x_0, x_1, x_2):

	# times to check
	times = [345.0, 714.0, 1079.0, 1446.0, 1935.0, 2018.0]
	sign = ""

	#data - time
	T = 2019.0			#total time
	num_steps = 2019	#number of steps
	dt = T/num_steps	#one step 

	#data 
	d_n = Constant(0.56)	#nutrients diffusion
	k = Constant(dt)		

	#fit no. 1
	e = 0.0764822198996
	M = 0.892319833114
	Ar = 1.87811473342
	Ag = 2.46062239811
	Ab = 2.14193159659

	#element, mesh
	print "Loading function space, geometry ..."
	nx = 20
	ny = 20
	nz = 20
	mesh = 	BoxMesh(Point(0,0,0), Point(6,6,6), nx, ny, nz)
	P1 = FiniteElement('P', tetrahedron, 1)
	element = MixedElement([P1, P1])
	V = FunctionSpace(mesh, element)
	
	#test function, solution
	v_1, v_2 = TestFunctions(V)
	u = Function(V)
	u_1, u_2 = split(u)

	#initial position of tumour
	sr = Point(x_0, x_1, x_2)
	
	#initial conditions 
	print "Loading initial conditions ..."

	class IniCond(Expression):

		def __init__(self, sr, a, **kwargs):
			self.sr = sr
			self.a = a

    		def eval(self, value, x):
				dx = (x[0] - self.sr[0])*(x[0] - self.sr[0])
				dy = (x[1] - self.sr[1])*(x[1] - self.sr[1])
				dz = (x[2] - self.sr[2])*(x[2] - self.sr[2])
				
				if (dx + dy + dz < 1.02):
					value[0] = a
				else:
					value[0] = 0.0
				
				value[1] = 1.17

		def value_shape(self):
        		return (2,)

	u_0 = IniCond(sr,a,degree=5)
	u_n = Function(V)
	u_n = interpolate(u_0, V)
	u_n1, u_n2 = split(u_n)

	#load blood vessels (our three types: red, green, blue)
	print "Loading blood vessels..."

	f = CellFunction("int", mesh, 3)

	if typ == "r":
		A = Ar
		with open("f_r.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
		
	elif typ == "g":
		A = Ag
		with open("f_g.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
	
	elif typ == "b":
		A = Ab
		with open("f_b.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
	
	r=0

	for cell in cells(mesh):
		f[r] = int(data[r])
		r = r+1

	class BloodVessels(Expression):

		def __init__(self, f, **kwargs):
			self.f = f
		
		def eval_cell(self, value, x, ufc_cell):
			value[0] = self.f[ufc_cell.index]

	f0 = BloodVessels(f, degree=5)

	#variational problem 
	F = ((u_1 - u_n1)*v_1/k)*dx + \
		d_c*dot(grad(u_1), grad(v_1))*dx - \
		g*u_1*u_2*v_1*dx  + \
		((u_2 - u_n2)*v_2/k)*dx + \
		d_n*dot(grad(u_2), grad(v_2))*dx + \
		d*u_2*v_2*u_1*dx + \
		e*u_2*v_2*(((u_2*u_2))/(M+((u_2*u_2))))*dx - \
		A*f0*v_2*dx

	#saving pictures
	vtk1 = File('cancer3D_c_'+str(x_0)+'_'+ str(x_1) + '_' + str(x_2) + '_' +str(typ) +  '/cell_density.pvd')
	vtk2 = File('cancer3D_s_'+str(x_0)+'_'+ str(x_1) + '_' + str(x_2) + '_' +str(typ) + '/food_density.pvd')
	
	#progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)
	
	#solution
	t = 0
	for n in range(num_steps):
		
		# update time
		t += dt
	
		# solution + update
		solve(F == 0, u)
	
		# save 
		_u_1, _u_2 = u.split()
		vtk1 << (_u_1, t)
		vtk2 << (_u_2, t)

		# compute volume 
		if t in times:
			prec = 0.1
			score = 0
			tresh = 16000.0
			for i in range(0, int((6/prec))):
				for j in range(0, int((6/prec))):
					for k in range(0, int((6/prec))):
						if _u_1(i*prec,j*prec,k*prec) >= tresh:
							score += 1 

			curr_vol = (score/(((6/prec))*((6/prec))*((6/prec))))*6*6*6
			sign = sign + str((curr_vol)) + str(" ")

		# update 
		u_n.assign(u)
		progress.update(t / T)

	return sign

def generate():

	a = 57436.0159738
	g = 0.0013034946847
	d = 1.59932679501e-08
	d_c = 9.95688081374e-05

	#type of vessels 
	typ_array = ["r", "g", "b"]
	typ_int = random.randrange(3)
	typ = typ_array[typ_int]
	
	#inital position
	x_0 = random.uniform(1.75,4.15)
	x_1 = random.uniform(1.75,4.15)
	x_2 = random.uniform(1.75,4.15)
	
	sign = rd_solver(a,g,d,d_c, typ, x_0, x_1, x_2)
	w = open("results.txt","a")
	w.write(str(a) + "\t" + str(g) + "\t" + str(d) + "\t" + str(d_c) + "\t" + str(x_0) + "\t" + str(x_1) + "\t" + str(x_2) + "\t" + str(sign) + "\n")
	w.close()	
	
generate()