import numpy as np
import scipy as sp


class BeamGaussian:
    def __init__(self, z, w0, n, lam, E0, lx, ly, Nx, Ny):

        """
        :param z: beam spreading distance
        :param w0: waist radius
        :param n: reflective index
        :param lam: wavelength
        :param E0: initial field
        :param lx: x-axis length of the screen
        :param ly: y-axis length of the screen
        :param Nx: x-axis discretization
        :param Ny: y-axis discretization
        """
        self.z: np.longdouble = z
        self.w0: np.longdouble = w0
        self.lam: np.longdouble = lam
        self.n: np.longdouble = n
        self.E0: np.longdouble = E0
        self.lx: np.longdouble = lx
        self.ly: np.longdouble = ly
        self.Nx: np.int_ = Nx
        self.Ny: np.int_ = Ny

    def field(self):
        k = 2 * np.pi * self.n / self.lam
        zr = np.pi * np.power(self.w0, 2) * self.n / self.lam   # Rayleigh range

        x = np.linspace(-self.lx/2, self.lx/2, self.Nx)
        y = np.linspace(-self.ly/2, self.ly/2, self.Ny)
        xv, yv = np.meshgrid(x, y)

        r2 = np.power(xv, 2) + np.power(yv, 2)

        w = self.w0 * np.sqrt(1 + np.power(self.z / zr, 2))    # hyperbolic relation
        R = self.z * (1 + np.power(zr / self.z, 2))    # radius of curvature

        psi = np.arctan(self.z / zr)    # Gouy phase

        field = self.E0 * (self.w0 / w) * np.exp(-r2 / np.power(w, 2)) * np.exp(-1j * (k * self.z + k * r2 / (2 * R) - psi))
        intensity = np.power(np.abs(field), 2)

        return intensity, field, xv, yv


class ThinLens:
    def __init__(self, D, field_input, f, lam, lx, ly):
        
        """
        :param D: lens diameter
        :param field_input: incident field
        :param f: lens focal length
        :param lam: wavelength
        :param lx: x-axis length of the screen
        :param ly: н-axis length of the screen
        """
             
        self.D: np.longdouble = D
        self.field_input: np.complex256 = field_input
        self.f: np.longdouble = f
        self.lam: np.longdouble = lam
        self.lx: np.longdouble = lx
        self.ly: np.longdouble = ly

    def output_field_forward(self):
        """
        Field after passing through a thin lens at forward beam propagation
        :return: output field, phase function for the lens, output intensity
        """

        k = 2*np.pi/self.lam

        Ny, Nx = np.shape(self.field_input)

        x = np.linspace(-self.lx/2, self.lx/2, Nx)
        y = np.linspace(-self.ly/2, self.ly/2, Ny)
        xv, yv = np.meshgrid(x, y)

        # transmission function(f>0 for collecting lens)
        lens = np.power(xv, 2)+np.power(yv, 2) <= np.power(self.D/2, 2)
        lens = lens.astype(float)
        phi = (-k/(2*self.f)*(np.power(xv, 2) + np.power(yv, 2)))*lens % (2*np.pi)
        phase = np.exp(1j*phi)

        field_output = self.field_input*phase
        intensity_output = np.power(np.abs(field_output), 2)

        return field_output, phi, intensity_output

    def output_field_backward(self):
        """
        Field after passing through a thin lens at backward beam propagation
        :return: output field, phase function for the lens, output intensity
        """

        k = 2 * np.pi / self.lam

        Ny, Nx = np.shape(self.field_input)

        x = np.linspace(-self.lx / 2, self.lx / 2, Nx)
        y = np.linspace(-self.ly / 2, self.ly / 2, Ny)
        xv, yv = np.meshgrid(x, y)

        # transmission function
        lens = np.power(xv, 2) + np.power(yv, 2) <= np.power(self.D / 2, 2)
        lens = lens.astype(float)
        phi = k / (2 * self.f) * (np.power(xv, 2) + np.power(yv, 2))*lens % (2 * np.pi)
        phase = np.exp(1j * phi)

        field_output = self.field_input * phase
        intensity_output = np.power(np.abs(field_output), 2)

        return field_output, phi, intensity_output


class CylindricalLens:
    def __init__(self, D, field_input, f, lam):

        """
        Field after passing through a cylindrical lens

        :param D: lens diameter
        :param field_input: incident field
        :param f: lens focal length
        :param lam: wavelength
        """

        self.D: np.longdouble = D
        self.field_input: np.complex256 = field_input
        self.f: np.longdouble = f
        self.lam: np.longdouble = lam

    def output_field(self):
        Ny = np.shape(self.field_input)[1]

        y = np.linspace(-self.D/2, self.D/2, Ny)

        yv = np.meshgrid(y, y)[1]
        # transmission function
        lens = np.absolute(yv) <= self.D/2
        lens = lens.astype(float)

        phase = np.exp(1j*(-2*np.pi/(2*self.lam*self.f)*np.power(yv, 2))*lens)

        field_output = self.field_input*lens*phase
        return field_output


class PropagationFresnel:
    def __init__(self, z, field_input, lam, lx, ly):
        """
        :param z: beam spreading distance
        :param field_input: input field
        :param lam: wavelength
        :param lx: x-axis length of the screen
        :param ly: y-axis length of the screen
        """

        self.z: np.longdouble = z
        self.field_input: np.complex256 = field_input
        self.lam: np.longdouble = lam
        self.lx: np.longdouble = lx
        self.ly: np.longdouble = ly

    def new_field(self):

        Ny, Nx = np.shape(self.field_input)

        x_nf = np.linspace(-self.lx/2, self.lx/2, Nx)
        y_nf = np.linspace(-self.ly/2, self.ly/2, Ny)

        # spatial frequencies
        kx_nf = sp.fft.fftfreq(Nx, np.diff(x_nf)[0]) * (2 * np.pi)
        ky_nf = sp.fft.fftfreq(Ny, np.diff(y_nf)[0]) * (2 * np.pi)
        kxv_nf, kyv_nf = np.meshgrid(kx_nf, ky_nf)

        k = 2*np.pi/self.lam

        field_input_fft = sp.fft.fft2(self.field_input, workers=8)

        # Fourier image of impulse response function
        ir_fft = -np.exp(1j*k*self.z)*np.exp(-1j*(self.lam*self.z/(4*np.pi))*(np.power(kxv_nf, 2)+np.power(kyv_nf, 2)))

        # Fourier image of output field
        field_output_fft = field_input_fft*ir_fft

        # output field
        field_output = sp.fft.ifft2(field_output_fft, workers=8)
        intensity_output = np.power(np.abs(field_output), 2)
        return field_output, intensity_output







