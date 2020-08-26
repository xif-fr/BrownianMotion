import numpy as np
import matplotlib.pyplot as plt
from math import *
import pysimul

def plot_means_synth (t, X, cmak=4):
	plt.figure(figsize=(15,5))
	plt.plot(t, X, lw=1, color='royalblue')
	if cmak != 0:
		plt.plot(t, pysimul.cma(X, cmak), color='navy')
	X_mean = np.mean(X)
	X_std = np.std(X)
	plt.axhline(X_mean, color='darkorange', label="Mean : {:.2e}".format(X_mean))
	plt.fill_between(t, X_mean-X_std, X_mean+X_std, facecolor='orange', alpha=0.3, label="Std. dev. : {:.1e}".format(X_std))
	(X_slope, X0) = np.polyfit(t, X, 1)
	plt.plot(t, X0 + X_slope*t, color='darkorange', lw=1, label="Change : {:.1e}".format(X_slope*(t[-1]-t[0])))
	plt.xlim([t[0],t[-1]])
	plt.legend(loc='upper right')
	return (X_mean, X_std, X_slope)

def check_gaussian (samples, bins, xlabel):
	plt.hist(samples, bins=50, density=True)
	x0 = np.mean(samples)
	sigma = np.std(samples)
	x = np.linspace(-5*sigma, +5*sigma, 100)
	plt.plot(x, np.exp(-((x-x0)/sigma)**2/2)/(sigma*np.sqrt(2*np.pi)), lw=2, linestyle='--', color='black')
	plt.xlabel(xlabel)
	return x0, sigma
