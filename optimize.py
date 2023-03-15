from math import sqrt, inf
import subprocess
import re
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


FILE_TO_OPTIMIZE = "MRSt_OMEGA"
PARAMETER_NAMES = ["Oct", "Q1", "Q2", "H1", "H2", "S1", "S2", "S3", "angle", "u1", "u2"]
MIN_OCT = -0.0005
MAX_OCT =  0.0005

try:
	with open(f"{FILE_TO_OPTIMIZE}_cache.pkl", "rb") as file:
		cache = pickle.load(file)
except FileNotFoundError:
	cache = {}


def system_quality(parameters):
	output = run_cosy(parameters)
	time_skew = float(re.search(r"Time skew \(ps/keV\) = +([-.\d]+)", output).group(1))
	tof_width = float(re.search(r"N 5 FPDESIGN Time Resol\.\(ps\) +([-.\d]+)", output).group(1))
	energy_width = float(re.search(r"N 4 FPDESIGN HO Resol\.RAY\(keV\) +([-.\d]+)", output).group(1))
	time_resolution = sqrt(tof_width**2 + (energy_width*time_skew)**2)
	print(f"this design has a time resolution of {time_resolution:.1f} ps and an energy resolution of {energy_width:.1f} keV")
	return 100*(time_resolution/100 + energy_width/400)


def get_defaults():
	with open(f"{FILE_TO_OPTIMIZE}.fox", "r") as f:
		script = f.read()
	values = []
	bounds = []
	for name in PARAMETER_NAMES[1:]:
		values.append(float(re.search(rf"{name} := ([-.\d]+);", script).group(1)))
		if name.startswith("S"):
			bounds.append((0, 2.0))
		elif name.startswith("angle"):
			bounds.append((0, 90))
		elif name.startswith("u"):
			bounds.append((-60, 60))
		else:
			bounds.append((-.05, .05))
	return np.array(values), bounds


def run_cosy(parameters):
	""" get the observable values at these perturbations """
	parameters = tuple(parameters)
	if parameters not in cache:
		with open(f'{FILE_TO_OPTIMIZE}.fox', 'r') as f:
			script = f.read()
		
		for i, name in enumerate(PARAMETER_NAMES):
			script = re.sub(rf"{name} := [-.\d]+;", f"{name} := {parameters[i]};", script)

		with open('temp.fox', 'w') as g:
			g.write(script)

		try:
			result = subprocess.run(['C:/Program Files/COSY 10.0/cosy.exe', 'temp'], capture_output=True, check=True)
		except subprocess.CalledProcessError:
			print(result.stdout.decode('ascii'))
			raise

		# store full parameter sets and their resulting COSY outputs in the cache
		output = result.stdout.decode('ascii')
		output = output[1036:]
		cache[parameters] = output

		with open(f"{FILE_TO_OPTIMIZE}_cache.pkl", "wb") as file:
			pickle.dump(cache, file)

	return cache[parameters]


def p_distance(*parameters):
	output = run_cosy(parameters)
	distance = float(re.search(r"N 3 FPDESIGN p-dist\(mm\) +([-.\d]+)", output).group(1))
	print(f"  {parameters[0]}T -> {distance:.2g}mm")
	return distance


def outer_objective(parameters):
	parameters = tuple(parameters)
	if parameters not in cache:
		try:
			Oct = optimize.brentq(p_distance, MIN_OCT, MAX_OCT, args=parameters, rtol=1e-5)
		except ValueError:
			print("You need to widen the octopole limits!")
			raise
			# Oct = min([MIN_OCT, MAX_OCT], key=lambda o: abs(p_distance(o, *parameters)))
		# store partial parameter sets and the corresponding best octopole strength in the cache
		cache[parameters] = Oct

	quality = system_quality((cache[parameters],) + parameters)
	print(f"{parameters} -> {quality:.2f}ps")
	return quality


def simplexify(x0, ranges):
	vertices = [np.array(x0)]
	for i in range(len(x0)):
		vertices.append(np.array(x0))
		vertices[-1][i] += (ranges[i][1] - ranges[i][0])/10


def optimize_design():
	defaults, bounds = get_defaults()
	initial_simplex = simplexify(defaults, bounds)
	result = optimize.minimize(
		outer_objective,
		defaults,
		bounds=bounds,
		method='Nelder-Mead',
		options=dict(initial_simplex=initial_simplex)
		)
	print(result.x)


if __name__ == '__main__':
	optimize_design()
