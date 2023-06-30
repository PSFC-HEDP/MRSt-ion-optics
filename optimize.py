import pickle
import re
import subprocess
from math import sqrt, inf, exp
from typing import Tuple, List, Union

import numpy as np
from numpy.typing import NDArray
from scipy import optimize


FILE_TO_OPTIMIZE = "MRSt_OMEGA"
# ORDER = 2
# PARAMETER_NAMES = ["Q2", "H2", "S1", "angle", "u1", "u2", "c1", "c2"]
ORDER = 3
PARAMETER_NAMES = ["Q2", "H2", "S1", "S2", "angle", "u1", "u2", "c1", "c2"]


with open(f'{FILE_TO_OPTIMIZE}.fox', 'r') as f:
	script = f.read()

try:
	with open(f"{FILE_TO_OPTIMIZE}_cache.pkl", "rb") as file:
		cache = pickle.load(file)
except FileNotFoundError:
	cache = {}


def optimize_design():
	""" optimize a COSY file by tweaking the given parameters to minimize the defined objective function """
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


def objective_function(parameters: List[float]) -> float:
	""" run COSY, read its output, and calculate a number that quantifies the system. smaller should be better """
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


def run_cosy(parameters: List[float]) -> str:
	""" get the observable values at these perturbations """
	parameters = tuple(parameters)
	if parameters not in cache or "### ERRORS IN CODE" in cache[parameters]:
		modified_script = script
		modified_script = re.sub(r"streamlined_mode := \d;", "streamlined_mode := 1;", modified_script)
		modified_script = re.sub(rf"order := \d;", f"order := {ORDER};", modified_script)
		for i, name in enumerate(PARAMETER_NAMES):
			modified_script = re.sub(rf"{name} := [-.\d]+;", f"{name} := {parameters[i]};", modified_script)

		with open('temp.fox', 'w') as f:
			f.write(modified_script)

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


def get_cosy_output(pattern: str, output: str) -> float:
	""" extract a single number from some COSY output using regex """
	match = re.search(pattern + r"([-.\de*]+)[\sa-z()/]*", output)
	if match is None:
		print(output)
		raise ValueError(f"couldn’t find /{pattern}/ in output")
	number = match.group(1)
	if "***" in number:
		return inf
	else:
		return float(number)


def get_defaults() -> Tuple[NDArray[float], List[Tuple[float, float]]]:
	""" read a COSY file to see what the free parameters are currently set to """
	values = []
	bounds = []
	for name in PARAMETER_NAMES:
		try:
			values.append(float(re.search(rf"{name} := ([-.\de]+);", script).group(1)))
		except AttributeError:
			print(script)
			raise RuntimeError(f"where is {name}?")
		if name.startswith("S"):
			bounds.append((0, 2.0))  # gaps must be positive, no more than 2m
		elif name.startswith("angle"):
			bounds.append((0, 90))  # don’t bend the beam more than 90°
		elif name.startswith("u"):
			bounds.append((-45, 45))  # entrance/exit plane angles should be within 45°
		elif name.startswith("c"):
			bounds.append((-5e-3, 5e-3))  # entrance/exit plane bending should be within 5mm
		elif name.startswith("H") or name.startswith("Q") or name.startswith("O"):
			bounds.append((-.05, .05))  # magnetic field strengths should be within 50mT
		else:
			raise ValueError(f"I don’t know what kind of quantity {name} is, so I can’t guess what its bounds should be.")
	return np.array(values), bounds


def simplexify(x0: Union[NDArray[float], List[float]],
               ranges: List[Tuple[float, float]]) -> NDArray[float]:
	""" build an initial simplex out of an initial guess, using their bounds as a guide """
	vertices = [np.array(x0)]
	for i in range(len(x0)):
		vertices.append(np.array(x0))
		vertices[-1][i] += (ranges[i][1] - ranges[i][0])/10
	return np.array(vertices)


if __name__ == '__main__':
	optimize_design()
