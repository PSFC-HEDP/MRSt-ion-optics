from math import sqrt
import subprocess
import re
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


time_skew = 8.1 # ps/keV
min_oct = -0.002
max_oct = -0.004

try:
	with open("cache.pkl", "rb") as file:
		cache = pickle.load(file)
except FileNotFoundError:
	cache = {}


def run_cosy(parameters):
	""" get the observable values at these perturbations """
	parameters = tuple(parameters)
	if parameters not in cache:

		cO, aQ, aH, bQ, bH = parameters

		with open('MRSt_final_70deg_proton.fox', 'r') as f:
			script = f.read()
		
		script = re.sub(r"AQ := [-.\d]+;", f"AQ := {aQ};", script)
		script = re.sub(r"AH := [-.\d]+;", f"AH := {aH};", script)
		script = re.sub(r"OCT := [-.\d]+;", f"OCT := {cO};", script)
		script = re.sub(r"BQ := [-.\d]+;", f"BQ := {bQ};", script)
		script = re.sub(r"BH := [-.\d]+;", f"BH := {bH};", script)

		with open('temp.fox', 'w') as g:
			g.write(script)

		try:
			result = subprocess.run(['C:/Program Files/COSY 10.0/cosy.exe', 'temp'], capture_output=True, check=True)
		except subprocess.CalledProcessError:
			print(result.stdout.decode('ascii'))
			raise

		cache[parameters] = result.stdout.decode('ascii')
		with open("cache.pkl", "wb") as file:
			pickle.dump(cache, file)

	return cache[parameters]


def p_distance(*parameters):
	output = run_cosy(parameters)
	distance = float(re.search(r"N 3 FPDESIGN p-dist\(mm\) +([-.\d]+)", output).group(1))
	print(f"  {parameters[0]} -> {distance:.2g}mm")
	return distance


def time_resolution(parameters):
	output = run_cosy(parameters)
	tof_width = float(re.search(r"N 5 FPDESIGN Time Resol\.\(ps\) +([-.\d]+)", output).group(1))
	energy_width = float(re.search(r"N 4 FPDESIGN HO Resol\.RAY\(keV\) +([-.\d]+)", output).group(1))
	total = sqrt(tof_width**2 + (energy_width*time_skew)**2)
	return total


def outer_objective(parameters):
	parameters = tuple(parameters)
	if parameters not in cache:
		try:
			cO = optimize.brentq(p_distance, min_oct, max_oct, args=parameters, rtol=1e-5)
		except ValueError:
			cO = min([min_oct, max_oct], key=lambda o: abs(p_distance(o, *parameters)))
		cache[parameters] = time_resolution((cO,) + parameters)

	print(f"{parameters} -> {cache[parameters]:.2f}ps")
	return cache[parameters]


def optimize_design():
	result = optimize.minimize(outer_objective, (-0.00485060696, -0.00052117182, -0.030908382, 0.00225), method='Nelder-Mead')
	print(result.x)


if __name__ == '__main__':
	optimize_design()
