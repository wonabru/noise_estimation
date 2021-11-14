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


if __name__ == '__main__':

	win_len = 400
	est = []
	for n in range(20):
		data = np.random.randn(win_len)
		est_noise, error = calc_noise_std(data)
		print('Noise estimation', est_noise, 'error:', error)
		est.append(est_noise)

	print('widnow length:', win_len, 'Avg. NTS:', np.mean(est), 'Error:', np.std(est))