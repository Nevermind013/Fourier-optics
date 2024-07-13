import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import major_classes as fo


import functions as func
import optimisation_algorithms as oa

# all sizes in mm
lx = 15.36*1e0
ly = 8.64*1e0
Lx = 15.36*1e0
Ly = 8.64*1e0
lam = 1064 * 1e-6

# discretisation
Nx = 1920
Ny = 1080

# z1 -- distance between source and DOE, z2 -- between DOE and lens, z3 -- between lens and screen
z1 = 400
z2 = 1000
z3 = 1000

T = 1.5*lam   # period of the grating
N = Lx / T

x0 = np.linspace(-lx/2, lx/2, Nx)
y0 = np.linspace(-ly/2, ly/2, Ny)
xv0, yv0 = np.meshgrid(x0, y0)

# just for plot
sfx = ScalarFormatter()
sfy = ScalarFormatter()
sfx.set_powerlimits((-int(lx/2*100), int(lx/2*100)))
sfy.set_powerlimits((-int(ly/2*100), int(ly/2*100)))


# transmission function for the 1d sin grating
DOE_surface = (np.abs(xv0) <= Lx/2)*(np.abs(yv0) <= Ly/2).astype(float)
# phase multiplier grating, created by DOE with 1d sin phase function
grating = np.exp(1j*np.power(np.sin(2*np.pi*xv0/T), 2))*DOE_surface

# initial field
E0 = (np.abs(xv0) <= Lx/20)*(np.abs(yv0) <= Ly/20).astype(float)

# from source to DOE
f1 = fo.PropagationFresnel(z=z1, field_input=E0, lam=lam, lx=lx, ly=ly)
E1, I_source = f1.new_field()

E_doe = E1*grating

# from DOE to lens
f2 = fo.PropagationFresnel(z=z2, field_input=E_doe, lam=lam, lx=lx, ly=ly)
E2 = f2.new_field()[0]

focus = z3
lens = fo.ThinLens(D=50, field_input=E2, f=focus, lam=lam, lx=lx, ly=ly)
E3 = lens.output_field_forward()[0]

# from lens to screen
f3 = fo.PropagationFresnel(z=z3, field_input=E3, lam=lam, lx=lx, ly=ly)
E4, intensity_screen = f3.new_field()

intensity_screen_1d = intensity_screen[int(Ny/2)]

# first-order selection mask
mask = (xv0 >= 3.4)*(xv0 <= 4.8)*(yv0 >= -0.85)*(yv0 <= 0.85).astype(float)

# define intensities for the optimising problem
intensity_target = intensity_screen*mask
intensity_source = np.power(np.abs(E1), 2)
intensity_target_1d = intensity_target[int(Ny/2)]

func.screen(xv=xv0, yv=yv0, zv=intensity_screen, title=r'$I_{screen}$ after 1d sin grating', sfx=sfx, sfy=sfy)
plt.show()
func.screen(xv=xv0, yv=yv0, zv=intensity_target, title=r'$I_{target}$', sfx=sfx, sfy=sfy)
plt.show()

# plot 1d intensity after propagating through the sin grating
fig, ax = plt.subplots(figsize=(7, 6), edgecolor='black', linewidth=3, frameon=True)
ax.plot(x0, intensity_screen_1d, color='blue', label=r'$I_{screen}^{1d}$')
ax.set_title(r'$I_{screen}^{1d}$')
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel('$intensity$')
ax.xaxis.set_major_formatter(sfx)
ax.yaxis.set_major_formatter(sfy)
ax.grid(which='both')
ax.legend()
plt.show()

# plot 1d target intensity
fig, ax = plt.subplots(figsize=(7, 6), edgecolor='black', linewidth=3, frameon=True)
ax.plot(x0, intensity_target_1d, color='blue', label=r'$I_{target}^{1d}$')
ax.set_title(r'$I_{target}^{1d}$')
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel('$intensity$')
ax.xaxis.set_major_formatter(sfx)
ax.yaxis.set_major_formatter(sfy)
ax.grid(which='both')
ax.legend()
plt.show()


phi_initial = np.random.rand(Ny, Nx)*2*np.pi

# phase_by_FD = oa.FastestDescent(intensity_target=intensity_target, intensity_source=intensity_source, z=z2, lx=lx,
#                                   ly=ly, lam=lam, accuracy=1e-3, phi_initial=phi_initial)
# phase_FD1d = phase_by_FD1d.phase_retrieval_1d()[0]
#
# func.screen(xv=xv0, yv=yv0, zv=phase_FD1d, title=r'$\phi(x,y)$', sfx=sfx, sfy=sfy)
# plt.show()

# phase_FD_2d = phase_by_FD.phase_retrieval_2d()[0]


phase_by_GS = oa.GerchbergSaxtonAlgorithm(intensity_target=intensity_target, intensity_source=intensity_source, z=z2,
                                          lx=lx, ly=ly, lam=lam, accuracy=1e-3, phi_initial=phi_initial)

phase_GS_1d = phase_by_GS.phase_retrieval_2d()[0]


# propagating from DOE to screen after optimise phase function
E_doe_new = E1*np.exp(1j*phase_GS_1d)

f2_new = fo.PropagationFresnel(z=z2, field_input=E_doe_new, lam=lam, lx=lx, ly=ly)
E2_new = f2_new.new_field()[0]

focus = z3
lens_new = fo.ThinLens(D=50, field_input=E2_new, f=focus, lam=lam, lx=lx, ly=ly)
E3_new = lens_new.output_field_forward()[0]

f3_new = fo.PropagationFresnel(z=z3, field_input=E3_new, lam=lam, lx=lx, ly=ly)
E4_new, intensity_new = f3_new.new_field()

delta = func.relative_error(intensity=intensity_new, target=intensity_target)
print('Relative error ='+str(delta))

func.screen(xv=xv0, yv=yv0, zv=intensity_target/np.sum(intensity_target), title=r'$I_{target}^{norm}$', sfx=sfx, sfy=sfy)
func.screen(xv=xv0, yv=yv0, zv=intensity_new/np.sum(intensity_new), title=r'$I_{numerical}^{norm}$', sfx=sfx, sfy=sfy)
plt.show()

fig, ax = plt.subplots(figsize=(8, 6), edgecolor='black', linewidth=3, frameon=True)
ax.plot(x0, intensity_new[int(Ny/2)]/np.sum(intensity_new[int(Ny/2)]), color='red', label=r'$I_{numerical}^{1d}$')
ax.plot(x0, intensity_target_1d/np.sum(intensity_target_1d), color='blue', label=r'$I_{target}^{1d}$')
ax.set_title(r'Results of one-dimensional optimisation')
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel('$intensity$')
ax.xaxis.set_major_formatter(sfx)
ax.yaxis.set_major_formatter(sfy)
ax.grid(which='both')
ax.legend()
plt.show()





