import numpy as np


def gaussian(x, a, b, c):
    return a * np.exp(-(x - b) ** 2 / (2 * c))


def parabola(x, a):
    return a * x ** 2


def func_g_norm(x, b):
    # Normalised Gaussian function
    # 1 Parameter
    return np.exp(-b * x**2)


def func_g_norm_e(b):
    return np.sqrt(1 / float(b))


def func_g_norm_e_err(b, b_error):
    return 0.5 * float(b)**(-1.5) * b_error


def func_e_norm(x, b):
    # Exponential function
    # 1 Parameter
    return np.exp(-b * x)


def func_e_norm_hwhm(b):
    return np.log(2) * 1 / float(b)


def func_e_norm_h_err(b, b_error):
    return np.log(2) * b**(-2) * b_error


def mean_func(freq, a, b):
	# return a * freq + b
	return a * np.power(freq, b)


def linear(x,a,b):
    return a * x + b


def linear_sym(x,a,b):
    return a * np.abs(x) + b

# def func5(x, para):
#     # Gaussian function
#     # 2 Parameter
#     return para[0]* np.exp(-para[1] * x**2)


# def func5_hwhm(para):
#     return np.sqrt(1 / float(para[1]))


# def func5_err(para, error):
#     return 0.5 * float(para[1])**(-1.5) * error[1]


# def func(x, a, b, c):
#     # Exponential function
#     # 3 Parameters
#     return a * np.exp(-b * x) + c


# def func3(x, a, b):
#     # Exponential function
#     # 2 Parameters
#     return a * np.exp(-b * x)
