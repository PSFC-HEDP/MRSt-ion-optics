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
UNITS = {'strength':'%', 'tiltx':'degrees', 'tilty':'degrees', 'shiftx':'cm', 'shifty':'cm', 'shiftz':'cm', 'length':'cm', 'aperture':'cm'}
OBSERVABLES = [ # the parameters we can measure (name, lower margin, upper margin, controller index)
	('Tilt Angle(deg)', -1, 1, 'strength_H2'),
	('p-dist(mm)', -1, 1, 'strength_O')]


def get_values(hexapole, octopole):
	""" get the observable values at these perturbations """
	with open('MRSt_tol.fox', 'r') as f:
		fox = f.read()
	for parameter in PARAMETERS:
		if parameter == 'strength_H2':
			fox = fox.replace(f'<<{parameter}>>', str(hexapole))
		elif parameter == 'strength_O':
			fox = fox.replace(f'<<{parameter}>>', str(octopole))
		if 'strength' in parameter:
			fox = fox.replace(f'<<{parameter}>>', '1')
		else:
			fox = fox.replace(f'<<{parameter}>>', '0')
	with open('temp.fox', 'w') as g:
		g.write(fox)
	try:
		res = subprocess.run(['cosy', 'temp'], capture_output=True, check=True)
	except subprocess.CalledProcessError:
		print(res.stdout.decode('ascii'))
	vals = []
	for line in res.stdout.decode('ascii').split('\r\n'):
		if 'FPDESIGN' in line:
			if "Tilt Angle" in line:
				vals.append(float(line.split()[-1]))
			elif "p-dist" in line:
				vals.append(float(line.split()[-1]))
	return vals


if __name__ == '__main__':
	x = np.linspace(0, 30, 12)
	y = np.linspace(-2, 2, 13)
	tilt = np.empty((y.size, x.size))
	bend = np.empty((y.size, x.size))
	for j in range(x.size):
		for i in range(y.size):
			print(i,j)
			tilt[i,j], bend[i,j] = get_values(x[j], y[i])
	plt.figure()
	plt.contourf(x, y, tilt)
	plt.title("tilt")
	plt.xlabel("hexapole strength")
	plt.ylabel("octopole strength")
	plt.colorbar()
	plt.figure()
	plt.contourf(x, y, bend)
	plt.title("bend")
	plt.xlabel("hexapole strength")
	plt.ylabel("octopole strength")
	plt.colorbar()
	plt.show()
