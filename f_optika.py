import numpy as np
import scipy as sp


class BeamGaussian():
    def __init__(self, z, w0, n, lam, E0, lx, ly, Nx, Ny):
        self.z: float = z
        # w0 -- spot size parameter
        self.w0: float = w0
        self.lam: float = lam
        self.n: float = n
        self.E0: float = E0
        self.lx: float = lx
        self.ly: float = ly
        self.Nx: int = Nx
        self.Ny: int = Ny

    def field(self):
        k = 2 * np.pi * self.n / self.lam
        zr = np.pi * np.power(self.w0, 2) * self.n / self.lam

        x = np.linspace(-self.lx/2, self.lx/2, self.Nx)
        y = np.linspace(-self.ly/2, self.ly/2, self.Ny)
        xv, yv = np.meshgrid(x, y)
        r2 = np.power(xv, 2) + np.power(yv, 2)

        w = self.w0 * np.sqrt(1 + np.power(self.z / zr, 2))
        R = self.z * (1 + np.power(zr / self.z, 2))

        psi = np.arctan(self.z / zr)

        E = self.E0 * (self.w0 / w) * np.exp(-r2 / np.power(w, 2)) * np.exp(-1j * (k * self.z + k * r2 / (2 * R) - psi))

        E_real = np.real(E)
        E_imag = np.imag(E)

        intensity = np.power(E_real, 2) + np.power(E_imag, 2)

        return intensity, E, xv, yv


class ThinLens():

    def __init__(self, D, E_input, f, lam, n, lx, ly):
        self.D: float = D
        self.E_input = E_input
        self.f: float = f
        self.lam: float = lam
        self.n: float = n
        self.lx: float = lx
        self.ly: float = ly

    def output_field(self):

        k = 2*np.pi*self.n/self.lam

        Nx: int = np.shape(self.E_input)[0]
        Ny: int = np.shape(self.E_input)[1]

        x = np.linspace(-self.lx/2, self.lx/2, Nx)
        y = np.linspace(-self.ly/2, self.ly/2, Ny)

        xv, yv = np.meshgrid(x, y)
        # функция пропускания
        lense = np.power(xv, 2)+np.power(yv, 2) <= np.power(self.D/2, 2)
        lense = lense.astype(float)

        phase = np.exp(-1j*(k/(2*self.f)*(np.power(xv, 2) + np.power(yv, 2))))*lense

        E_output = self.E_input*phase*lense

        return E_output


class CylindricalLens():
    def __init__(self, D, E_input, f, lam, n):
        self.D: float = D
        self.E_input: np.ndarray = E_input
        self.f: float = f
        self.lam: float = lam
        self.n: float = n

    def output_field(self):
        Ny: int = np.shape(self.E_input)[1]

        y = np.linspace(-self.D/2, self.D/2, Ny)

        yv = np.meshgrid(y, y)[1]
        # функция пропускания
        lense = np.absolute(yv) <= self.D/2
        lense = lense.astype(float)

        phase = np.exp(1j*(-2*np.pi*self.n/(2*self.lam*self.f)*np.power(yv, 2))*lense)
        E_output = self.E_input*lense*phase
        return E_output


class PropogationFresnel():
    def __init__(self, z, field, lam, lx, ly):
        self.z: float = z
        self.field: np.ndarray = field
        self.lam: float = lam
        self.lx: float = lx
        self.ly: float = ly

    def new_field(self):
        Nx: int = np.shape(self.field)[0]
        Ny: int = np.shape(self.field)[1]

        x_nf: np.ndarray = np.linspace(-self.lx/2, self.lx/2, Nx)
        y_nf: np.ndarray = np.linspace(-self.ly/2, self.ly/2, Ny)
        xv_nf, yv_nf = np.meshgrid(x_nf, y_nf)

        k = 2*np.pi/self.lam

        # фурье-преобразование входящего поля
        E_input_fft = sp.fft.fft2(self.field)

        # пространственные частоты
        kx_nf = sp.fft.fftfreq(Nx, np.diff(x_nf)[0])*(2*np.pi)
        ky_nf = sp.fft.fftfreq(Ny, np.diff(y_nf)[0])*(2*np.pi)
        kxv_nf, kyv_nf = np.meshgrid(kx_nf, ky_nf)

        # фурье-образ импульсного ответа
        # ir_fft = np.exp(1j*self.z*2*np.pi/self.lam)*np.exp(-1j*np.pi*self.lam*self.z*(np.power(kxv_nf, 2)+np.power(kyv_nf, 2)))
        ir_fft = np.exp(-1j*self.z*np.sqrt(np.power(k, 2)-np.power(kxv_nf,2)-np.power(kyv_nf, 2)))

        # фурье-образ поля в плоскости z=const
        E_output_fft = E_input_fft*ir_fft

        # поле в плоскости z=const
        E_output = sp.fft.ifft2(E_output_fft)

        E_output_real = np.real(E_output)
        E_output_imag = np.imag(E_output)

        intensity_output = np.power(E_output_real, 2)+np.power(E_output_imag, 2)

        return intensity_output, xv_nf, yv_nf, E_input_fft, kxv_nf, kyv_nf, E_output_fft, E_output







