import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import major_classes as fo
import GB_algorithm as gb
from scipy import interpolate

params = {'font.family': 'Times New Roman', 'legend.fontsize': 18, 'figure.figsize': (9.5, 6), 'axes.labelsize': 18, 'axes.titlesize': 20,
          'xtick.labelsize': 18, 'ytick.labelsize': 18, 'figure.titlesize': 20}
plt.rcParams.update(params)
plt.style.use(['science', 'notebook'])
sf = ScalarFormatter()



path_gb = r'C:\Users\gunne\Fourier-optics\GB\gaussian beam'
path_out = r'C:\Users\gunne\Fourier-optics\GB\output'
path_delta = r'C:\Users\gunne\Fourier-optics\GB\delta'
def screen(xv, yv, mass, title, sfx,sfy, name):

    fig, ax = plt.subplots(figsize=(7, 7), edgecolor='black', linewidth=3, frameon=True)
    im = ax.pcolormesh(xv, yv, mass, cmap='inferno')
    ax.set_title(title)
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')

    ax.xaxis.set_major_formatter(sfx)
    ax.yaxis.set_major_formatter(sfy)
    fig.colorbar(im)
    # plt.savefig(name+'.png')
    # ax.text(0., 0., 'BAM', backgroundcolor='red', label='1111')
    return None


def optical_system(lx, ly, Nx, Ny, lam, w0, I_target,phi_initial, z0, z1):
    ''' Создаём расчётную сетку'''
    x0 = np.linspace(-lx / 2, lx / 2, 2500)
    y0 = np.linspace(-ly / 2, ly / 2, 2500)
    xv0, yv0 = np.meshgrid(x0, y0)


    ''' Создаем источник'''
    beam1 = fo.BeamGaussian(z=5*1e-2, w0=w0, n=1, lam=lam, E0=1, lx=lx, ly=ly, Nx=2500, Ny=2500)
    I_gb = beam1.field()[0]
    E_gb = beam1.field()[1]
    phase0 = (np.angle(E_gb) + 2*np.pi) % (2*np.pi)
    intensity_source = np.copy(I_gb)
    screen(xv0, yv0, intensity_source, r'Профиль интенсивности в исходной плоскости', sf, sf, '1')

    f1 = fo.PropogationFresnel(z=1000*1e-2, field=E_gb, lam=lam, lx=lx, ly=ly)
    I = f1.new_field()[0]
    screen(xv0, yv0, I, r'Профиль интенсивности на экране', sf, sf,'1')
    plt.show()













    # screen(xv0, yv0, I_target, r'$I_{target}$', sf, '1')
    # plt.show()

    ''' ГС для поиска фазовой функции ДОЭ'''
    phase1 = gb.GB_algorithm(I_target, intensity_source, z1, lx=lx, ly=ly, lam=lam, accuracy=1e-3,
                             initial=phi_initial, M=256)

    phi1, intensity_target_new = phase1.phase()
    # phi1 = phi1 - phase0

    ''' Поле после прохождения ДОЭ'''
    E_output = E_gb * np.exp(1j * phi1)

    f1 = fo.PropogationFresnel(z=z1, field=E_output, lam=lam, lx=lx, ly=ly)
    E1 = f1.new_field()[7]

    focus = z1
    doe0 = fo.ThinLens(4 * 1e-2, E1, focus, lam, lx, ly)
    E = doe0.output_field()[0]

    f2 = fo.PropogationFresnel(z=z1, field=E, lam=lam, lx=lx, ly=ly)
    intensity_out = f2.new_field()[0]
    E_out = f2.new_field()[7]
    delta_target = np.sum(I_target)
    delta_f = np.sqrt(np.sum(np.power(np.absolute(E_out) - np.sqrt(I_target), 2)) / delta_target)


    # screen(xv0, yv0, intensity_out, r'$I~after~GS$', sf, 1)
    # screen(xv0, yv0, phi1, r'$Phase function$', sf, 1)
    # plt.show()

    return delta_f, phi1






if __name__ == "__main__":
    ''' Параметры системы'''
    lx = 15.36*1e-3
    ly = 15.36*1e-3
    lam = 1064*1e-9
    w0 = 1*1e-3
    z0 = 15*1e-2
    z1 = 10*1e-2
    l = 8.64*1e-3
    # Дискретизация
    Nx = 1920
    Ny = 1080
    # Масштабирование экрана
    sf.set_powerlimits((-lx/2*100, lx/2*100))

    x0 = np.linspace(-lx/2, lx/2, Nx)
    y0 = np.linspace(-ly/2, ly/2, Ny)
    xv0, yv0 = np.meshgrid(x0, y0)
    # искомое распределение
    u0 = (np.abs(xv0) < l / 100) * (np.abs(yv0) < l / 10)
    u0 = u0.astype(float) * 4
    u01 = (np.abs(yv0) < l / 100) * (np.abs(xv0) < l / 10)
    u01 = u01.astype(float) * 4
    u01 = u0 + u01

    # delta_f = optical_system(lx=lx, ly=ly, Nx=Nx, Ny=Ny, lam=lam, w0=w0, I_target=u01, phi_initial=2*np.pi*np.random.rand(Ny, Nx), z0=z0, z1=z1)
    # print('Относительная ошибка:'+str(delta_f))


    dimx = np.array([320, 520, 640, 720, 1280, 1600, 1920])
    dimy = np.array([180, 292, 360, 405, 720, 900, 1080])
    N = 7
    error = np.array([])
    count= 0

    phi_init = 2 * np.pi * np.random.rand(dimy[0], dimx[0])
    for n in range(N-1):

        Nx = dimx[n]
        Ny = dimy[n]
        x0 = np.linspace(-lx / 2, lx / 2, Nx)
        y0 = np.linspace(-ly / 2, ly / 2, Ny)
        xv0, yv0 = np.meshgrid(x0, y0)
        # искомое распределение
        u0 = (np.abs(xv0) < l / 100) * (np.abs(yv0) < l / 10)
        u0 = u0.astype(float) * 4
        u01 = (np.abs(yv0) < l / 100) * (np.abs(xv0) < l / 10)
        u01 = u01.astype(float) * 4
        u01 = u0 + u01
        delta_f, phi1 = optical_system(lx=lx, ly=ly, Nx=Nx, Ny=Ny, lam=lam, w0=w0, I_target=u01,
                                 phi_initial=phi_init, z0=z0, z1=z1)
        error = np.append(error, delta_f)
        func = interpolate.interp2d(x0, y0, phi1, kind='cubic')

        x = np.linspace(-lx / 2, lx / 2, dimx[n+1])
        y = np.linspace(-ly / 2, ly / 2, dimy[n + 1])
        phi_init = np.copy(func(x, y))
        count += 1
        print(count)

    # dim0 = np.linspace(100, 1500, 100)
    # dim0 = [int(n) for n in dim0]

    # error0 = np.array([0.49442859141164586,
    #                    0.4263098266844485,
    #                    0.37426493786620396,
    #                    0.37303003473846025,
    #                    0.3945618182377807,
    #                    0.35568484057320265,
    #                    0.3494186569923503,
    #                    0.3339135243781279,
    #                    0.37488218239035737,
    #                    0.38439234906861786,
    #                    0.3552881312487875,
    #                    0.3832205989894711,
    #                    0.38477505397139156,
    #                    0.35788001183372287,
    #                    0.3491916574848808,
    #                    0.3728696407804173,
    #                    0.3526645182721876,
    #                    0.3587044383630889,
    #                    0.37178007808337554,
    #                    0.35509228242097196,
    #                    0.3583146623729653,
    #                    0.35435276592411175,
    #                    0.3772358545027789,
    #                    0.37137192941936126,
    #                    0.3801409440552141,
    #                    0.37041867577733667,
    #                    0.36503612585404577,
    #                    0.36538647341710656,
    #                    0.3630596480272804,
    #                    0.37350785296898315,
    #                    0.36346700012389804,
    #                    0.3760318170757677,
    #                    0.3855283089188848,
    #                    0.37954696171754626,
    #                    0.3801993961455561,
    #                    0.3710200726943328,
    #                    0.3909241897268212,
    #                    0.38374979740842996,
    #                    0.3873255294607186,
    #                    0.3747499721421436,
    #                    0.3789106398162288,
    #                    0.37107331990575676,
    #                    0.3645994139891241,
    #                    0.3766130455782126,
    #                    0.36281736165517847,
    #                    0.3849079739667293,
    #                    0.38299413273396216,
    #                    0.3947762006807487,
    #                    0.37905070416593556,
    #                    0.3932567232775726,
    #                    0.39341864178546115,
    #                    0.3724640524163783,
    #                    0.38336099418200126,
    #                    0.3884975820233043,
    #                    0.37407407778807333,
    #                    0.3930785509340505,
    #                    0.40107793479953313,
    #                    0.39021238879566256,
    #                    0.3886870144490888,
    #                    0.3966720904451646,
    #                    0.39091588259396404,
    #                    0.3925430032365225,
    #                    0.39984542732624145,
    #                    0.3968664814851046,
    #                    0.39906292662176196,
    #                    0.41765061166546763,
    #                    0.42292081436095197,
    #                    0.42489483529212,
    #                    0.42401982508025854,
    #                    0.4343057780914188,
    #                    0.4259846185924541,
    #                    0.42481736476481274,
    #                    0.4165886146030543,
    #                    0.41315790329907226,
    #                    0.4195440538175817,
    #                    0.4276652763685989,
    #                    0.4205428854164272,
    #                    0.42132206701462294,
    #                    0.41598026565100493,
    #                    0.41369282647568717,
    #                    0.41863063636157494,
    #                    0.409132350256339,
    #                    0.4105848688934324,
    #                    0.4182113046720529,
    #                    0.40652981783399944,
    #                    0.4076321725757645,
    #                    0.40300724483352324,
    #                    0.40171624820526614,
    #                    0.40444809058342923,
    #                    0.39570477152226285,
    #                    0.4064546437370084,
    #                    0.4038658462089886,
    #                    0.39196150023997894,
    #                    0.3973093194287653,
    #                    0.40323539745352593,
    #                    0.40593241411242975,
    #                    0.4015374178435561,
    #                    0.4172419351446914,
    #                    0.4087355814331284,
    #                    0.4074102019628095
    #                    ])
    fig, ax = plt.subplots(figsize=(7, 6), edgecolor='black', linewidth=3, frameon=True)
    ax.set_title(r'$Error(N)$')
    ax.set_xlabel(r'$N_x$')
    ax.set_ylabel('$Error$')
    # ax.scatter(dim0, error0, s=50, color='blue', label='random')
    ax.scatter(dimx[:N-1], error, s=70, color='red', label='cubic')
    ax.legend()
    plt.show()
    # plt.savefig(name+'.png')
    # ax.text(0., 0., 'BAM', backgroundcolor='red', label='1111')













    ''' Распространение с преобразованием фазовой функции'''
    # from scipy import interpolate
    #
    # func = interpolate.interp2d(x0, y0, phi1, kind='linear')
    #
    #
    #
    #
    #
    # Nx = 600
    # Ny = 600
    #
    # x0 = np.linspace(-l / 2, l / 2, Nx)
    # y0 = np.linspace(-l / 2, l / 2, Ny)
    # xv0, yv0 = np.meshgrid(x0, y0)
    #
    # phi_in = func(x0, y0)
    #
    # # щель размерами D*d
    # a = l
    # u0 = (np.abs(xv0) < l / 100) * (np.abs(yv0) < l / 10)
    # u0 = u0.astype(float) * 4
    #
    # u01 = (np.abs(yv0) < l / 100) * (np.abs(xv0) < l / 10)
    # u01 = u01.astype(float) * 4
    #
    # u01 = u0 + u01
    #
    # # screen(xv0, yv0, u0, r'rectangular slit, $a=$'+str(a)+'m', sf, 1)
    #
    # # круглое отверстие диаметром D
    # u1 = (np.power(xv0, 2) + np.power(yv0, 2) <= np.power(l / 5, 2)) * (
    #             np.power(xv0, 2) + np.power(yv0, 2) >= np.power(l / 7, 2))
    # u1 = u1.astype(float) * 4
    # # screen(xv0, yv0, u1, 'Circular hole', sf, 1)
    #
    # u2 = (np.power(xv0, 2) + np.power(yv0, 2) <= np.power(0.3 * 1e-3, 2))
    # u2 = u2.astype(float) * 0.01
    # # screen(xv0, yv0, u2, 'Circular hole', sf, 1)
    #
    # z_gb = 20 * 1e-2
    # z_gb1 = 200 * 1e-2
    # z = 30 * 1e-2
    #
    # # # пучок Гаусса
    #
    # beam1 = fo.BeamGaussian(z=z, w0=1 * 1e-3, n=1, lam=lam, E0=1, lx=l, ly=l, Nx=Nx, Ny=Ny)
    # I_gb = beam1.field()[0]
    # E_gb = beam1.field()[1]
    #
    # xv_gb = beam1.field()[2]
    # yv_gb = beam1.field()[3]
    #
    # phase0 = np.angle(E_gb)
    #
    # doe_func = (np.power(xv0, 2) + np.power(yv0, 2) <= np.power(l / 9, 2))
    # doe_func = doe_func.astype(float)
    #
    # beam2 = fo.BeamGaussian(z=0.01 * z, w0=1 * 1e-3, n=1, lam=lam, E0=1, lx=l, ly=l, Nx=Nx, Ny=Ny)
    # I_gb2 = beam2.field()[0]
    #
    # # интенсивность источника и интенсивность на экране
    # intensity_source = np.copy(I_gb)
    # intensity_target = np.copy(I_gb2)
    #
    # screen(xv_gb, yv_gb, intensity_source, r'$I_{source}$', sf, '1')
    # screen(xv_gb, yv_gb, intensity_target, r'$I_{target}$', sf, '1')
    # # screen(xv_gb, yv_gb, doe_func, r'$DOE$', sf, '1')
    # plt.show()
    #
    # phase1 = gb.GB_algorithm(intensity_target, intensity_source, z, lx=l, ly=l, lam=lam, accuracy=1e-3,
    #                          initial=phi_in, M=256)
    #
    # phi1, intensity_target_new = phase1.phase()
    # phi1 = phi1 - phase0
    #
    # E_output = E_gb * np.exp(1j * phi1)
    #
    # f1 = fo.PropogationFresnel(z=z, field=E_output, lam=lam, lx=l, ly=l)
    # E1 = f1.new_field()[7]
    #
    # focus = z
    # doe0 = fo.ThinLens(4 * 1e-2, E1, focus, lam, l, l)
    # E = doe0.output_field()[0]
    #
    # f2 = fo.PropogationFresnel(z=z, field=E, lam=lam, lx=l, ly=l)
    # intensity_out = f2.new_field()[0]
    # E_out = f2.new_field()[7]
    #
    # delta_target = np.sum(intensity_target)
    #
    # delta_f = np.sqrt(np.sum(np.power(np.absolute(E_out) - np.sqrt(intensity_target), 2)) / delta_target)
    #
    # print('Intensity error:' + str(delta_f))
    #
    # screen(xv0, yv0, intensity_out, r'$I~after~GS~2$', sf, 1)
    # screen(xv0, yv0, phi1, r'$PF$', sf, 1)
    # plt.show()


























