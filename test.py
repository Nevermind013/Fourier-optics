import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import f_optika as fo
params = {'font.family': 'Times New Roman', 'legend.fontsize': 18, 'figure.figsize': (9.5, 6), 'axes.labelsize': 18, 'axes.titlesize': 20,
          'xtick.labelsize': 18, 'ytick.labelsize': 18, 'figure.titlesize': 20}
plt.rcParams.update(params)
plt.style.use(['science', 'notebook'])
sf = ScalarFormatter()


def fft(kxv, kyv, mass, title):

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.pcolormesh(sp.fft.fftshift(kxv), sp.fft.fftshift(kyv), np.abs(sp.fft.fftshift(mass)))
    ax.set_title(title)
    ax.set_xlabel('$k_x$ [mm$^{-1}$]')
    ax.set_ylabel('$k_y$ [mm$^{-1}$]')

    # plt.xlim(-1e5, 1 * 1e5)
    # plt.ylim(-1e5, 1*1e5)
    plt.xlim(-np.amax(kxv)*0.1, np.amax(kxv)*0.1)
    plt.ylim(-np.amax(kyv)*0.1, np.amax(kyv)*0.1)
    plt.show()
    return None

def screen(xv, yv, mass, title, sf):

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pcolormesh(xv, yv, mass, cmap='inferno')
    ax.set_title(title)
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')

    ax.xaxis.set_major_formatter(sf)
    ax.yaxis.set_major_formatter(sf)
    plt.show()
    return None


if __name__ == "__main__":
    # параметры(m)
    D = 0.05*1e-3
    d = 5*1e-3
    # длина волны(m)
    lam = 632 * 1e-9
    # размер экрана(m)
    l = 4*1e-2
    # Дискретизация
    Nx = 7000
    Ny = 7000

    # Масштабирование экрана
    sf.set_powerlimits((-l/2*100, l/2*100))

    x0 = np.linspace(-l/2, l/2, Nx)
    y0 = np.linspace(-l/2, l/2, Ny)
    xv0, yv0 = np.meshgrid(x0, y0)

    # щель размерами D*d
    u0 = (np.abs(xv0) < D/2) * (np.abs(yv0) < d)
    u0 = u0.astype(float)
    # screen(xv0, yv0, u0, 'rectangular slit', sf)


    # круглое отверстие диаметром 3D
    u1 = np.power(xv0, 2)+np.power(yv0, 2) <= np.power(3*D/2, 2)
    u1 = u1.astype(float)
    screen(xv0, yv0, u1, 'Circular hole', sf)

    # пучок Гаусса
    z_gb = 10*1e-2
    beam1 = fo.BeamGaussian(z=z_gb, w0=0.7 * 1e-3, n=1, lam=lam, E0=1, lx=l, ly=l, Nx=Nx, Ny=Ny)

    I_gb = beam1.field()[0]

    E_gb = beam1.field()[1]

    xv_gb = beam1.field()[2]
    yv_gb = beam1.field()[3]

    # screen(xv_gb, yv_gb, I_gb, r'Gaussian Beam before lens, $z=$'+str(z_gb)+' m', sf)

    f1 = fo.PropogationFresnel(z=10 * 1e-2, field=u1, lam=lam, lx=l, ly=l)

    intensity = f1.new_field()[0]
    E_output = f1.new_field()[7]
    # E_input_fft = f1.new_field()[3]
    # E_output_fft = f1.new_field()[6]
    # kxv = f1.new_field()[4]
    # kyv = f1.new_field()[5]

    # fft(kxv, kyv, E_input_fft, r'FFT $E_{input}$')
    # fft(kxv, kyv, E_output_fft, r'FFT $E_{output}$')
    # интенсивность в плоскости перед линзой
    screen(xv0, yv0, intensity, r'$I_{output}$', sf)

    # график интенсивности в плоскости перед линзой
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d')
    axes.set_title(r'$I_{output}(x,y)$')
    axes.plot_surface(xv0, yv0, intensity)
    plt.show()

    doe1 = fo.ThinLens(4 * 1e-2, E_output, 10 * 1e-2, lam, 2, l, l)
    # doe2 = fo.CylindricalLens(4*10-2, E_output, 30*1e-2, lam, 2)

    # поле сразу после линзы
    E1 = doe1.output_field()

    f2 = fo.PropogationFresnel(z=10 * 1e-2, field=E1, lam=lam, lx=l, ly=l)

    intensity2 = f2.new_field()[0]
    E_output2 = f2.new_field()[7]
    E_input_fft2 = f2.new_field()[3]
    E_output_fft2 = f2.new_field()[6]
    kxv2 = f2.new_field()[4]
    kyv2 = f2.new_field()[5]
    # картина на экране
    screen(xv0, yv0, intensity2, r'$I_{output}~on~a~screen$', sf)

    # График интенсивности в плоскости экрана
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d')
    axes.set_title(r'$I_{output}(x,y)~on~a~screen$')
    axes.plot_surface(xv0, yv0, intensity2)
    # plt.show()

