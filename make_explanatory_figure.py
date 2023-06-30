from math import pi, sin, cos, ceil, nan

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from numpy.typing import NDArray

r_wire = .1

def main():
	fig, axs = plt.subplots(3, 2, figsize=(4.2, 6), facecolor="none")
	for i, ax in enumerate(axs[:, 0]):
		visible_xs = []
		visible_ys = []
		present_xs = []
		present_ys = []
		sines = []
		poles = 2*(i + 1)
		width = 6/poles
		radius = 1
		if width > 2:
			radius /= width/2
			width = 2
		thickness = ceil(width/2/(2*r_wire))*(2*r_wire)
		# θ = np.linspace(0, 2*pi, 144, endpoint=False)
		# ax.plot(radius*np.cos(θ), radius*np.sin(θ), "k--", linewidth=1.2)
		patches = []
		for j in range(poles):
			angle = (j + 0.5)*2*pi/poles
			sine = j%2*2 - 1
			ax.fill(*rotate(angle,
			                np.array([radius, radius, radius + thickness, radius + thickness]),
			                np.array([-width/2, width/2, width/2, -width/2])),
			        color="#999", zorder=2)
			ax.plot(*rotate(angle, np.array([radius, radius]), np.array([-width/2, width/2])),
			        color="k", linewidth=1.2, zorder=3)
			for side in [-1, 1]:
				num_wires_to_show = round(thickness/(2*r_wire))
				for k in range(4*num_wires_to_show):
					x = radius + (k + 0.5)*2*r_wire
					y = sine*side*width/2
					x, y = rotate(angle, x, y)
					if k < num_wires_to_show:
						patches.append(Circle((x, y), r_wire, facecolor="white", edgecolor="k", linewidth=1.2))
						if side < 0:
							patches.append(Circle((x, y), .4*r_wire, facecolor="k", edgecolor="none"))
						else:
							d = .4*r_wire
							ax.plot([x - d, x + d], [y - d, y + d], "k-", linewidth=0.8, zorder=4)
							ax.plot([x - d, x + d], [y + d, y - d], "k-", linewidth=0.8, zorder=4)
						visible_xs.append(x)
						visible_ys.append(y)
					present_xs.append(x)
					present_ys.append(y)
					sines.append(side)
		ax.add_collection(PatchCollection(patches, match_original=True, zorder=3))
		r_image = np.max(np.abs(visible_xs + visible_ys)) + r_wire
		r_aperture = np.min(np.hypot(visible_xs, visible_ys))

		X, Y = np.meshgrid(np.linspace(-r_aperture, r_aperture, 101),
		                   np.linspace(-r_aperture, r_aperture, 101), indexing="ij")
		A = np.zeros_like(X)
		for x0, y0, sine in zip(present_xs, present_ys, sines):
			R = np.hypot(X - x0, Y - y0)
			A += sine/R
		A[np.hypot(X, Y) > r_aperture] = nan
		ax.contour(X, Y, A, levels=np.linspace(-20, 20, 16), colors="#700", linewidths=.6,
		           linestyles="solid", negative_linestyles="solid", zorder=1)

		ax.set_xlim(-r_image, r_image)
		ax.set_ylim(-r_image, r_image)
		ax.set_aspect('equal', 'box')
		ax.set_axis_off()

	for ax in axs[:, 1]:
		ax.set_axis_off()
	plt.tight_layout()
	plt.show()

def rotate(angle: float, x: NDArray[float] | float, y: NDArray[float] | float
           ) -> tuple[NDArray[float] | float, NDArray[float] | float]:
	return (cos(angle)*x - sin(angle)*y,
	        sin(angle)*x + cos(angle)*y)

if __name__ == "__main__":
	main()
