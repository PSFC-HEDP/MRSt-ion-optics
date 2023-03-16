from math import sqrt, inf
import subprocess
import re
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


FILE_TO_OPTIMIZE = "MRSt_OMEGA"
PARAMETER_NAMES = ["Q1", "Q2", "H1", "H2", "S1", "S2", "S3", "angle", "u1", "u2"]
# FILE_TO_OPTIMIZE = "MRSt_OMEGA_quadratic"
# PARAMETER_NAMES = ["Q1", "Q2", "H1", "H2", "S1", "S2", "angle", "u1", "u2"]
# FILE_TO_OPTIMIZE = "MRSt_OMEGA_linear"
# PARAMETER_NAMES = ["Q1", "Q2", "S1", "S2", "angle", "u1", "u2"]

with open(f'{FILE_TO_OPTIMIZE}.fox', 'r') as f:
	script = f.read()

try:
	with open(f"{FILE_TO_OPTIMIZE}_cache.pkl", "rb") as file:
		cache = pickle.load(file)
except FileNotFoundError:
	cache = {}


def optimize_design():
	defaults, bounds = get_defaults()
	initial_simplex = simplexify(defaults, bounds)
	result = optimize.minimize(
		system_quality,
		defaults,
		bounds=bounds,
		method='Nelder-Mead',
		options=dict(initial_simplex=initial_simplex)
		)
	print(result)


def system_quality(parameters):
	output = run_cosy(parameters)
	time_skew = float(re.search(r"Time skew \(ps/keV\) += +([-.\d]+)", output).group(1))
	tof_width = float(re.search(r"FPDESIGN Time Resol\.\(ps\) +([-.\d]+)", output).group(1))
	energy_width = float(re.search(r"FPDESIGN HO Resol\.RAY\(keV\) +([-.\d]+)", output).group(1))
	time_resolution = sqrt(tof_width**2 + (energy_width*time_skew)**2)
	print(f"this design has a time resolution of {time_resolution:.1f} ps and an energy resolution of {energy_width:.1f} keV")

	quality = 100*(time_resolution/100 + energy_width/400)
	print(f"{parameters} -> {quality:.2f}ps")
	return quality


def run_cosy(parameters):
	""" get the observable values at these perturbations """
	parameters = tuple(parameters)
	if parameters not in cache or "### ERRORS IN CODE" in cache[parameters]:
		for i, name in enumerate(PARAMETER_NAMES):
			modified_script = re.sub(rf"{name} := [-.\d]+;", f"{name} := {parameters[i]};", script)

		with open('temp.fox', 'w') as g:
			g.write(modified_script)

		try:
			result = subprocess.run(['C:/Program Files/COSY 10.0/cosy.exe', 'temp'], capture_output=True, check=True)
		except subprocess.CalledProcessError:
			print(result.stdout.decode('ascii'))
			raise

		# store full parameter sets and their resulting COSY outputs in the cache
		output = result.stdout.decode('ascii')
		if "### ERROR" in output:
			print(re.sub(r"[\n\r]+", "\n", output))
			raise RuntimeError("COSY threw an error")
		output = output[1036:]
		cache[parameters] = output

		with open(f"{FILE_TO_OPTIMIZE}_cache.pkl", "wb") as file:
			pickle.dump(cache, file)

	return cache[parameters]


def get_defaults():
	values = []
	bounds = []
	for name in PARAMETER_NAMES:
		values.append(float(re.search(rf"{name} := ([-.\d]+);", script).group(1)))
		if name.startswith("S"):
			bounds.append((0, 2.0))
		elif name.startswith("angle"):
			bounds.append((0, 90))
		elif name.startswith("u"):
			bounds.append((-45, 45))
		else:
			bounds.append((-.05, .05))
	return np.array(values), bounds


def simplexify(x0, ranges):
	vertices = [np.array(x0)]
	for i in range(len(x0)):
		vertices.append(np.array(x0))
		vertices[-1][i] += (ranges[i][1] - ranges[i][0])/10


if __name__ == '__main__':
	optimize_design()
