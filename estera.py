import numpy as np


def entropy(data: np.array, eps: float) -> float:
	multiplier = (1.0 - 0.564) * eps
	list_theta = []
	for i in range(len(data)):
		for l1 in range(1, len(data) - i):
			theta = 0
			p = np.nan
			for p in range(len(data) - i - l1):
				theta += abs(data[i + p] - data[i + l1 + p])
				if theta >= multiplier * (p + 1):
					break

			if p <= 0:
				p = np.nan
			if i + l1 + p < len(data):
				list_theta.append(p)

	mean_theta = np.nanmean(list_theta)
	if mean_theta > 1:
		k2 = np.log(mean_theta / (mean_theta - 1))
	else:
		k2 = 0
	return k2


def k2_noisy_calc(data: np.array) -> tuple:
	eps_max = np.std(data)
	eps = np.arange(0, eps_max, eps_max / 100)
	return eps, np.array(list(map(lambda x: entropy(data, x), eps)))


def calc_noise_std(data: np.array):

	_std = np.std(data)
	_mean = np.mean(data)

	data = (data - _mean) / _std * 0.7

	eps, k2 = k2_noisy_calc(data)

	def fun(sigma):
		p = 0.3441717 - 1 / np.log(sigma)
		no_eps = np.argmin(np.abs(eps - sigma))
		no_sigma = np.argmax(eps ** p * k2)
		# print(abs(no_sigma - no_eps))
		return abs(no_sigma - no_eps)

	min_ret = 100
	est_sig = np.nan
	for sig in np.arange(0.001, 1, 0.001):
		ret = fun(sig)
		if min_ret > ret:
			min_ret = ret
			est_sig = sig
		if min_ret == 0:
			break

	return est_sig / 0.7 * _std, min_ret


def henon_attractor(x, y, a=1.4, b=0.3):
	'''
	Computes the next step in the Henon
	map for arguments x, y with kwargs a and
	b as constants.
	'''
	x_next = 1 - a * x ** 2 + y
	y_next = b * x
	return x_next, y_next


if __name__ == '__main__':

	win_len = 1000
	est = []
	level_of_noise = 0.1
	for n in range(10):
		data = np.random.randn(win_len)
		# starting point
		X = []
		Y = []
		x_next = 0
		y_next = 0
		# add points to array
		for i in range(win_len):
			x_next, y_next = henon_attractor(x_next, y_next)
			X.append(x_next)
			Y.append(y_next)
		# Normalize X
		X = (np.array(X) - np.mean(X)) / np.std(X)
		data = (np.array(data) - np.mean(data)) / np.std(data)

		# add noise
		X = X + level_of_noise * data

		# estimate noise
		est_noise, error = calc_noise_std(X)
		print('Theoretical noise level:', np.std(level_of_noise * data), 'Noise estimated from data', est_noise, 'error:', error)
		est.append(est_noise)

	print('Theoretical noise level:', level_of_noise, 'widnow length:', win_len, 'Avg. estimated NTS:', np.mean(est), 'Error:', np.std(est))