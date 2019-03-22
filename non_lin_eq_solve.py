import numpy as np

def eq_solve(f,fx,fy,fz,g,gx,gy,gz,h,hx,hy,hz,x0,y0,z0,err = 1e-5,maxiter = 50):

	x1 = x0
	y1 = y0
	z1 = z0
	
	for i in range(0, maxiter):
	  A = np.array([[fx(x1,y1,z1), fy(x1,y1,z1), fz(x1,y1,z1)],
		     [gx(x1,y1,z1),gy(x1,y1,z1),gz(x1,y1,z1)],
		     [hx(x1,y1,z1),hy(x1,y1,z1),hz(x1,y1,z1)]])

  	  B = np.array([[-f(x1,y1,z1)],[-g(x1,y1,z1)],[-h(x1,y1,z1)]])
	  C = np.linalg.solve(A,B)
	  if abs(f(x1,y1,z1)) <= err and abs(g(x1,y1,z1)) <= err and abs(h(x1,y1,z1))<= err:
		
	  	return (x1,y1,z1)
	  else:
	  	x1 = x0 + float(C[0])
        	y1 = y0 + float(C[1])
	       	z1 = z0 + float(C[2])
		x0 = x1
		y0 = y1
		z0 = z1

	if i == (maxiter-1):
		print "Could not find solution in", i," iterations"

