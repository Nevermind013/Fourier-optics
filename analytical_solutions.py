import numpy as np
import scipy as sp


class SquareFresnel:
    def __init__(self, z, lam, lx, ly, Nx, Ny, a):
        """
        analytical solution for a screen with a square hole in the Fresnel approximation

        :param z: beam spreading distance
        :param lam: wavelength
        :param lx: x-axis length of the screen
        :param ly: н-axis length of the screen
        :param Nx: x-axis discretization
        :param Ny: y-axis discretization
        :param a: square side
        """
        self.z: float = z
        self.lam: float = lam
        self.lx: float = lx
        self.ly: float = ly
        self.Nx: int = Nx
        self.Ny: int = Ny
        self.a: float = a

    def intensity(self):

        k = 2*np.pi/self.lam
        x = np.linspace(-self.lx / 2, self.lx / 2, self.Nx)
        y = np.linspace(-self.ly / 2, self.ly / 2, self.Ny)
        xv, yv = np.meshgrid(x, y)

        intensity_theory = np.zeros((self.Nx, self.Ny))

        def feta1(x):
            return -np.sqrt(k / (np.pi * self.z)) * (self.a / 2 + x)

        def feta2(x):
            return np.sqrt(k / (np.pi * self.z)) * (self.a / 2 - x)

        def func_cos(t):
            return np.cos(np.pi * np.power(t, 2) / 2)

        def func_sin(t):
            return np.sin(np.pi * np.power(t, 2) / 2)

        def fc(v):
            return sp.integrate.quad(func_cos, 0, v)[0]

        def fs(y):
            return sp.integrate.quad(func_sin, 0, y)[0]

        for i in range(self.Nx):
            for j in range(self.Ny):
                ksi1 = feta1(xv[i, j])
                ksi2 = feta2(xv[i, j])

                eta1 = feta1(yv[i, j])
                eta2 = feta2(yv[i, j])

                intensity_theory[i, j] = ((1/4)*(np.power((fc(ksi2) - fc(ksi1)), 2)+np.power((fs(ksi2) - fs(ksi1)), 2))
                                          * (np.power((fc(eta2) - fc(eta1)), 2)+np.power((fs(eta2) - fs(eta1)), 2)))

        return intensity_theory


class CircleFraunhofer:
    def __init__(self, z, lam, lx, ly, Nx, Ny, D):
        """
        Analytical solution for a screen with a circular hole in the Fraunhofer approximation (something wrong, to finish)

        :param z: beam spreading distance
        :param lam: wavelength
        :param lx: x-axis length of the screen
        :param ly: y-axis length of the screen
        :param Nx: x-axis discretization
        :param Ny: y-axis discretization
        :param D: circle diameter
        """
        self.z: float = z
        self.lam: float = lam
        self.lx: float = lx
        self.ly: float = ly
        self.Nx: int = Nx
        self.Ny: int = Ny
        self.D: float = D

    def intensity(self):
        k = 2 * np.pi / self.lam
        x = np.linspace(-self.lx / 2, self.lx / 2, self.Nx)
        y = np.linspace(-self.ly / 2, self.ly / 2, self.Ny)
        xv, yv = np.meshgrid(x, y)

        intensity_theory = np.zeros((self.Nx, self.Ny))

        for i in range(self.Nx):
            for j in range(self.Ny):
                r0 = np.sqrt(np.power(xv[i, j], 2) + np.power(yv[i, j], 2))
                intensity_theory[i, j] = np.power(k*self.D*self.D/(8*self.z), 2)*np.power((4*self.z*sp.special.jv(1, (
                        k*self.D*r0)/(2*self.z))/(k*self.D*r0)), 2)

        return intensity_theory











