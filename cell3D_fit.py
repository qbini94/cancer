#biblioteki
import numpy as np
from fenics import * 
import time
from pyswarm import pso

def rd_solver(a, g, d, d_c, typ):

	#cost function and times to check
	times = {345.0:6.83, 714.0:13.6, 1079.0:19.13, 1446.0:26.86}
	cost = 0
	sign = ""

	#data - time
	T = 2019.0			#czas symulacji
	num_steps = 2019	#ilosc krokow
	dt = T/num_steps	#jeden przedzial

	#data to simulation
	d_n = Constant(0.56)	#dyfuzja dla skladnikow odzywcz
	k = Constant(dt)		#sie przyda

	#parameters from fit no. 1
	e = 0.0764822198996
	M = 0.892319833114
	Ar = 1.87811473342
	Ag = 2.46062239811
	Ab = 2.14193159659

	#mesh, element ...
	print "Loading function space, geometry ..."
	nx = 20
	ny = 20
	nz = 20
	mesh = 	BoxMesh(Point(0,0,0), Point(6,6,6), nx, ny, nz)
	P1 = FiniteElement('P', tetrahedron, 1)
	element = MixedElement([P1, P1])
	V = FunctionSpace(mesh, element)
	
	#test functions, solutions
	v_1, v_2 = TestFunctions(V)
	u = Function(V)
	u_1, u_2 = split(u)

	#initial position of cancer
	sr = Point(3.0,2.5,4.1)

	#warunki poczatkowe, nowa wersja
	print "Loading initial conditions ..."

	class IniCond(Expression):

		def __init__(self, sr, a, **kwargs):
			self.sr = sr
			self.a = a

    		def eval(self, value, x):
				dx = (x[0] - self.sr[0])*(x[0] - self.sr[0])
				dy = (x[1] - self.sr[1])*(x[1] - self.sr[1])
				dz = (x[2] - self.sr[2])*(x[2] - self.sr[2])
				
				#initial radius and concentration
				if (dx + dy + dz < 1.02):
					value[0] = a
				else:
					value[0] = 0.0
				
				#initial food_density
				value[1] = 1.17

		def value_shape(self):
        		return (2,)

	u_0 = IniCond(sr,a,degree=5)
	u_n = Function(V)
	u_n = interpolate(u_0, V)
	u_n1, u_n2 = split(u_n)

	#load blood vessels
	print "Loading blood vessels..."

	f = CellFunction("int", mesh, 3)

	if typ == "r":
		A = Ar
		with open("f_r.txt") as p:
			data =  (p.read().splitlines())[0]
		p.close()
		
	else:
		print "Zly typ naczyn!"
		
	#here are different types of vessels used in example
	"""
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
	"""
	
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

	#variational formulation 
	F = ((u_1 - u_n1)*v_1/k)*dx + \
		d_c*dot(grad(u_1), grad(v_1))*dx - \
		g*u_1*u_2*v_1*dx  + \
		((u_2 - u_n2)*v_2/k)*dx + \
		d_n*dot(grad(u_2), grad(v_2))*dx + \
		d*u_2*v_2*u_1*dx + \
		e*u_2*v_2*(((u_2*u_2))/(M+((u_2*u_2))))*dx - \
		A*f0*v_2*dx

	#saving pictures (disabled)
	#vtk1 = File('cancer3D_c_'+str(a)+'_'+ str(g) + '_' + str(d) + '_' +str(d_c) + '/cell_density.pvd')
	#vtk2 = File('cancer3D_s_'+str(a)+'_'+ str(g) + '_' + str(d) + '_' +str(d_c) + '/food_density.pvd')
	
	#progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)
	
	#solution 
	t = 0
	for n in range(num_steps):
		
		# update time
		t += dt
	
		# solve + update
		solve(F == 0, u)
	
		# save and save pictures (disabled)
		_u_1, _u_2 = u.split()
		#vtk1 << (_u_1, t)
		#vtk2 << (_u_2, t)

		# compute cost 
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
			cost += (times[t] - curr_vol)*(times[t] - curr_vol)
			sign = sign + str((times[t] - curr_vol)) + str(" ")
			
		# update
		u_n.assign(u)
		progress.update(t / T)

	return cost, sign


def cost_function_red(x):

	#parameters to fit
	a = x[0]
	g = x[1]
	d = x[2]
	d_c = x[3]

	cost, sign = rd_solver(a,g,d,d_c,"r")
	w = open("results.txt","a")
	w.write(str(a) + "\t" + str(g) + "\t" + str(d) + "\t" + str(d_c) + "\t" + str(cost) + "\t" + str(sign) + "\n")
	w.close()
	return cost

xopt, fopt = pso(cost_function_red, [54000.0,0.00115,0.000000008,0.00009],[70000.0,0.0015,0.00000005,0.00015], debug = True)