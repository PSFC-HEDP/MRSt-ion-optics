import pickle
import re
import subprocess
from math import sqrt, inf, exp

import numpy as np
from scipy import optimize

FILE_TO_OPTIMIZE = "MRSt_OMEGA_quadratic"
PARAMETER_NAMES = ["Q1", "Q2", "H1", "H2", "S1", "S2", "angle", "u1", "u2"]
# FILE_TO_OPTIMIZE = "MRSt_OMEGA"
# PARAMETER_NAMES = ["Q1", "Q2", "H1", "H2", "S1", "S2", "S3", "angle", "u1", "u2"]

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
		objective_function,
		defaults,
		bounds=bounds,
		method='Nelder-Mead',
		options=dict(initial_simplex=initial_simplex)
		)
	print(result)


def objective_function(parameters):
	output = run_cosy(parameters)
	time_skew = get_cosy_output(r"Time skew \(ps/keV\) += +", output)
	tof_width = get_cosy_output(r"FPDESIGN Time Resol\.\(ps\) +", output)
	energy_width = get_cosy_output(r"FPDESIGN HO Resol\.RAY\(keV\) +", output)
	tilt_angle = get_cosy_output(r"FPDESIGN Tilt Angle\(deg\) +", output)
	time_resolution = sqrt(tof_width**2 + (energy_width*time_skew)**2)

	cost = 150*(time_resolution/150 + energy_width/300 + (exp(abs(tilt_angle) - 80)/5))
	print("[", end="")
	for parameter in parameters:
		print(f"{parameter:.6g},", end="")
	print(f"]\n\t->   {time_resolution:.2f}ps + {energy_width:.2f}keV + {tilt_angle:.2f}° = {cost:.2f}ps")
	return cost


def run_cosy(parameters):
	""" get the observable values at these perturbations """
	parameters = tuple(parameters)
	if parameters not in cache or "### ERRORS IN CODE" in cache[parameters]:
		modified_script = re.sub(r"streamlined_mode := \d;", "streamlined_mode := 1;", script)
		for i, name in enumerate(PARAMETER_NAMES):
			modified_script = re.sub(rf"{name} := [-.\d]+;", f"{name} := {parameters[i]};", modified_script)

		with open('temp.fox', 'w') as g:
			g.write(modified_script)

		result = subprocess.run(['C:/Program Files/COSY 10.0/cosy.exe', 'temp'], capture_output=True, check=True)
		if result.returncode > 0:
			print(result.stdout.decode('ascii'))
			raise RuntimeError()

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


def get_cosy_output(pattern, output):
	match = re.search(pattern + r"([-.\d*]+)", output)
	if match is None:
		print(output)
		raise ValueError(f"couldn’t find /{pattern}/ in output")
	number = match.group(1)
	if "***" in number:
		return inf
	else:
		return float(number)


def get_defaults():
	values = []
	bounds = []
	for name in PARAMETER_NAMES:
		try:
			values.append(float(re.search(rf"{name} := ([-.\de]+);", script).group(1)))
		except AttributeError:
			print(script)
			raise RuntimeError(f"where is {name}?")
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
	return np.array(vertices)


if __name__ == '__main__':
	optimize_design()
