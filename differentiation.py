def richardson(f, df, x, h):
	der = (1/(12.0*h)) *(f(x - 2*h) - 8*f(x - h) + 8*f(x + h) -  f(x + 2*h))

	true = df(x)
	error = abs(true - der)

	return (der,error)
