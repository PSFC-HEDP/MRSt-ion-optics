import numpy as np
import matplotlib.pyplot as plt
import subprocess

MAKE_GRAFS = False

FINNESS = 0.1**(3/3) # the amount we care about the tolerances being exact
assert FINNESS < 1
PARAMETERS = [
		'shiftz_H1', 'shiftz_O', 'shiftz_H2', 'shiftx_H1', 'shiftx_O', 'shiftx_H2', 'shifty_H1', 'shifty_O', 'shifty_H2',
		'length_H1', 'length_O', 'length_H2',
		'tiltx_H1', 'tiltx_O', 'tiltx_H2', 'tilty_H1', 'tilty_O', 'tilty_H2',
		'strength_Q1', 'strength_H1', 'strength_O', 'strength_Q2', 'strength_H2',
		'aperture_distance'] # the parameters we can change
OBSERVABLES = [ # the parameters we can measure (name, lower margin, upper margin, controller index)
	('HO Resol.RAY(keV)', -np.inf, .2, None),
	('Time Resol.(ps)', -np.inf, 3, None),
	('y-Size(mm)', -np.inf, .8, None),
	('Plane Length(m)', -.03, .03, None),
	('Tilt Angle(deg)', -1, 1, 'strength_H2'),
	('p-dist(mm)', -1, 1, 'strength_O')]

UNITS = {'strength':'%', 'tiltx':'degrees', 'tilty':'degrees', 'shiftx':'cm', 'shifty':'cm', 'shiftz':'cm', 'length':'cm', 'aperture':'cm'}


def basis_vec(idx, val):
	""" return a vector whare the only nonzero element is v[idx] = val """
	v = np.zeros(len(PARAMETERS))
	v[idx] = val
	return v


def get_values(X):
	""" get the observable values at these perturbations """
	X = np.array(X)
	assert X.size == len(PARAMETERS), f"{X.size} != {len(PARAMETERS)}"
	with open('MRSt_tol.fox', 'r') as f:
		fox = f.read()
	for parameter, x in zip(PARAMETERS, X):
		unit = get_units(parameter)
		if unit == '%':
			formatted = f'{1 + x/100:f}'
		elif unit == 'cm':
			formatted = f'{x/100:f}'
		else:
			formatted = f'{x:f}'
		fox = fox.replace(f'<<{parameter}>>', formatted)
	with open('temp.fox', 'w') as g:
		g.write(fox)
	try:
		res = subprocess.run(['cosy', 'temp'], capture_output=True, check=True)
	except subprocess.CalledProcessError:
		print(res.stdout.decode('ascii'))
	vals = np.full(len(OBSERVABLES), np.nan)
	for line in res.stdout.decode('ascii').split('\r\n'):
		if 'FPDESIGN' in line:
			for j, (key, lo, hi, controller) in enumerate(OBSERVABLES):
				if key in line:
					vals[j] = float(line.split()[-1])
					break
	return vals


def is_acceptable(y, y_min, y_max, i=None):
	""" does this set of parameters push one of the observable
		engineering features out of bounds?
		if `i != None`, then don't worry about observables that are
		being controlled by parameters other than `i`.
	"""
	for j in range(y_min.size): # check each observable
		if OBSERVABLES[j][3] is None or OBSERVABLES[j][3] == PARAMETERS[i]: # that is not being corrected by a different parameter
			if y[j] < y_min[j] or y[j] > y_max[j]: # to see if it is out of bounds
				return False # if so, it's not acceptable
	return True # if we make it this far, we're good


def find_tolerance(initial_guess, y_min, y_max, i):
	""" find the largest number less than `initial_guess` that gets the
		`i`th observable in bounds. the sampled perturbations will all
		be of the same sign as `initial_guess`.
		returns: the positive tolerance
	"""
	tol = initial_guess
	while True: # until I say you're done,
		x = basis_vec(i, tol)
		y = get_values(x) # try out this change
		if is_acceptable(y, y_min, y_max, i): # if everything is in bounds
			return abs(tol) # take it and go
		else: # if something is wrong
			tol *= 0.1**(1/3) # reduce the tolerance and try again


def get_name(param_code):
	return param_code.replace('_', ' of ').replace('Q', 'quadrupole ').replace('H', 'hexapole ').replace('O', 'the octupole')

def get_units(param_code):
	return UNITS[param_code.split('_')[0]]


if __name__ == '__main__':
	x0 = np.zeros(len(PARAMETERS)) # get the base parameters
	y0 = get_values(x0) # get the base observables by running the program
	y_min = np.array([y + lo_bound for y, (_, lo_bound, _, _) in zip(y0, OBSERVABLES)])
	y_max = np.array([y + hi_bound for y, (_, _, hi_bound, _) in zip(y0, OBSERVABLES)])

	tol = np.array(3*[.5] + 6*[.02] + 3*[.1] + 6*[.1] + 5*[.5] + [1]) # begin finding the tolerances
	slopes = np.full((x0.size, y0.size), np.nan) # and the direction of the dependencies
	for i in range(x0.size):
		# tol_plus = find_tolerance(10, y_min, y_max, i) # look up
		# tol_minus = find_tolerance(-tol_plus, y_min, y_max, i) # look down
		# tol[i] = min(tol_plus, tol_minus)*FINNESS # take the more restrictive one

		ξ = np.array([-2*tol[i], -tol[i], 0, tol[i], 2*tol[i]]) # now plot out the dependency
		υ = np.array([get_values(basis_vec(i, perturbation)) for perturbation in ξ])
		slopes[i,:] = (υ[3,:] - υ[1,:])/(ξ[3] - ξ[1])

		for j in range(y0.size):
			plt.figure()
			# if np.isfinite(OBSERVABLES[j][1]):
			# 	plt.plot(ξ, [y0[j] + OBSERVABLES[j][1]]*len(ξ), 'w-')
			# if np.isfinite(OBSERVABLES[j][2]):
			# 	plt.plot(ξ, [y0[j] + OBSERVABLES[j][2]]*len(ξ), 'w-')
			if OBSERVABLES[j][3] is None or OBSERVABLES[j][3] == PARAMETERS[i]:
				plt.plot(ξ, υ[:,j], f'C{j}-')
			else:
				plt.plot(ξ, υ[:,j], f'C{j}--')
			plt.xlabel(f"{PARAMETERS[i]} ({get_units(PARAMETERS[i])})")
			plt.ylabel(OBSERVABLES[j][0])
			plt.savefig(f"figs/dependency-{i}-{j}.png")
		# plt.show()
		plt.close('all')

		print(f"The tolerance on the {get_name(PARAMETERS[i])} is ±{tol[i]:.1g} {get_units(PARAMETERS[i])}")

	print()

	tol_for_each_type_of_param = dict()
	for i in range(x0.size):
		type_of_param = PARAMETERS[i].split('_')[0]
		tol_for_each_type_of_param[type_of_param] = min(tol_for_each_type_of_param.get(type_of_param, np.inf), tol[i])
	for i in range(x0.size):
		type_of_param = PARAMETERS[i].split('_')[0]
		tol[i] = tol_for_each_type_of_param[type_of_param]

	print()

	for j in range(y0.size): # now check the compound effects
		x = x0 - tol*np.sign(slopes[:,j])
		minimum = get_values(x)[j]
		x = x0 + tol*np.sign(slopes[:,j])
		maximum = get_values(x)[j]
		print(f"The {OBSERVABLES[j][0]} range is {minimum} < {y0[j]} < {maximum}")
