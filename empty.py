import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
sf = ScalarFormatter()
#
# # амплитуда
# E0 = 1
# n = 1
lam = 660 * 1e-9
# w0 = 2 * lam
#
# k = 2 * np.pi * n / lam
# zr = np.pi * np.power(w0, 2) * n / lam
#
# N = 1000
#
# x = np.linspace(-70 * lam, 70 * lam, N)
# y = np.linspace(-70 * lam, 70 * lam, N)
# xv, yv = np.meshgrid(x, y)
# r2 = np.power(xv, 2) + np.power(yv, 2)
#
# r = np.ones((N, N))
# for i in range(N):
#     for j in range(N):
#         if (xv[i, j] > 0 and yv[i, j] > 0) or (xv[i, j] < 0 and yv[i, j] < 0):
#             r[i, j] = np.sqrt(xv[i, j]**2+yv[i, j]**2)
#         else:
#             r[i, j] = -np.sqrt(xv[i, j] ** 2 + yv[i, j] ** 2)
#
#
# z = np.linspace(-100*lam, 100*lam, N)
# def GB(z):
#     w = w0 * np.sqrt(1 + np.power(z / zr, 2))
#     R = z * (1 + np.power(zr / z, 2))
#
#     psi = np.arctan(z / zr)
#
#
#     E = E0*(w0/w)*np.exp(-r2/np.power(w, 2))*np.exp(-1j*(k*z+k*r2/(2*R)-psi))
#
#
#     I = np.real(E)**2+np.imag(E)**2
#     return I
#
#
# I_mass = np.array([])
#
# # for i in range(N):
# #     I_mass = np.append(I_mass, GB(z[i]))
#
# print(np.shape(I_mass))
#
#
# print('zr='+str(zr))
# z =10*1e-2
# w = w0 * np.sqrt(1 + np.power(z / zr, 2))
# print('w='+str(w))
#
# I = GB(50*lam)
# fig, ax = plt.subplots(figsize=(7, 6))
# im = ax.pcolormesh(xv/lam, yv/lam, I)
# fig.colorbar(im, ax=ax, label=r'$I/I_0$')
# ax.set_title('I1')
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$')
# plt.show()
#
#
# # fig = plt.figure()
# # axes = fig.add_subplot(projection='3d')
# # axes.set_title(r'$I$')
# # # axes.plot_surface(sp.fft.fftshift(kxv), sp.fft.fftshift(kyv), np.abs(sp.fft.fftshift(E_ft)))
# # axes.plot_surface(zv1/lam, r1/lam, I1)
# # plt.show()

num = -np.pi*1.5/(lam*10*1e-2)
print(num)