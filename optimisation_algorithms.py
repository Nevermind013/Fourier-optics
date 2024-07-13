import numpy as np
import scipy as sp
import autograd.numpy as np
from autograd import elementwise_grad as egrad
import torch

import major_classes as fo
import functions as func


class FastestDescent:
    def __init__(self, intensity_target, intensity_source, z, lx, ly, lam, accuracy, phi_initial, levels=256):
        """
                :param intensity_target: target intensity profile
                :param intensity_source: intensity profile in the DOE plane
                :param z: beam spreading distance between DOE and screen
                :param lx: x-axis length of the screen
                :param ly: y-axis length of the screen
                :param lam: wavelength
                :param accuracy: optimization termination criterion
                :param phi_initial: initial approximation of the DOE phase function
                :param levels: number of quantization levels
        """
        self.intensity_target: np.ndarray = intensity_target
        self.intensity_source: np.ndarray = intensity_source
        self.z: float = z
        self.lx: float = lx
        self.ly: float = ly
        self.lam: float = lam
        self.accuracy: float = accuracy
        self.phi_initial: np.ndarray = phi_initial
        self.levels: int = levels

    def phase_retrieval_1d(self):
        """
        Optimisation algorithm for 1d phase grating

        :return: phase function and binary mask for DOE
        """
        k = 2 * np.pi / self.lam
        Ny, Nx = np.shape(self.intensity_target)

        amplitude_source = np.sqrt(self.intensity_source)[int(Ny / 2)]
        intensity_target_1d = self.intensity_target[int(Ny / 2)]

        x_nf: np.ndarray = np.linspace(-self.lx / 2, self.lx / 2, Nx)
        y_nf: np.ndarray = np.linspace(-self.ly / 2, self.ly / 2, Ny)

        y = y_nf[int(Ny / 2)]

        region = (intensity_target_1d > np.amax(intensity_target_1d) / 5).astype(float)

        dx = self.lx / Nx
        dy = self.ly / Ny

        # array of spatial frequencies
        kx_nf = sp.fft.fftfreq(Nx, dx) * (2 * np.pi)
        ky_nf = sp.fft.fftfreq(Ny, dy) * (2 * np.pi)

        ky = ky_nf[int(Ny / 2)]

        # impulse response function for direct beam propagation
        impulse_response_straight = -np.exp(1j * k * self.z) * np.exp(
            -1j * (self.z / (2 * k)) * (np.power(kx_nf, 2) + np.power(ky, 2)))

        focus = self.z
        # transmission function
        lens = np.power(x_nf, 2) + np.power(y, 2) <= np.power(50 / 2, 2)
        lens = lens.astype(float)

        # phase function for thin lens with focal distance named focus
        phi_lens = (-k / (2 * focus) * (np.power(x_nf, 2) + np.power(y, 2))) % (2 * np.pi)
        phase = np.exp(1j * phi_lens * lens)

        target_norm = np.sum(intensity_target_1d)
        target_normalised = intensity_target_1d/target_norm

        def target_function(phase_func):
            field_n = np.fft.ifft((np.fft.fft(
                np.fft.ifft(np.fft.fft(amplitude_source*np.exp(1j*phase_func)) * impulse_response_straight) * phase) *
                impulse_response_straight))
            intensity = np.power(np.abs(field_n), 2)
            intensity_norm = np.sum(intensity)
            intensity_normalised = intensity/intensity_norm

            return np.sum(np.power(target_normalised*region - intensity_normalised*region, 2))

        phi_n = np.copy(self.phi_initial)[int(Ny/2)]
        eps_n = target_function(phi_n)
        count = 0   # iteration counter
        # iter_arr = np.array([])
        # eps_arr = np.array([])

        dcdphi = egrad(target_function)

        dcdphi_n = dcdphi(phi_n)    # gradient at the point phi_n

        s_n = -eps_n/np.power(np.linalg.norm(dcdphi_n, ord=2), 2)

        while True:
            phi_new = phi_n + s_n*dcdphi_n

            if np.abs(eps_n) <= self.accuracy:
                break
            else:
                # iter_arr = np.append(iter_arr, count)
                count += 1
                # eps_arr = np.append(eps_arr, eps_n)
                print("iteration:" + str(count))
                phi_n = np.copy(phi_new)
                eps_n = target_function(phi_new)
                dcdphi_n = dcdphi(phi_new)

                s_n = -eps_n / np.power(np.linalg.norm(dcdphi_n, ord=2), 2)

        # plt.loglog(iter_arr, eps_arr)
        # plt.show()
        phi_new = phi_new % (2*np.pi)

        step = 2*np.pi/self.levels

        phi_1d = (phi_new//step)*step
        phi_binary_1d = (phi_new//step)

        phi = np.array([])
        phi_binary = np.array([])

        for i in range(Ny):
            phi = np.append(phi, phi_1d)
            phi_binary = np.append(phi_binary, phi_binary_1d)

        phi = np.reshape(phi, (Ny, Nx))
        phi_binary = np.reshape(phi_binary, (Ny, Nx))

        return phi, phi_binary

    def phase_retrieval_2d(self):
        """
        Optimisation algorithm for 2d phase grating

        :return: 2d phase function and 2d binary mask for DOE
        """

        region = (self.intensity_target > np.amax(self.intensity_target) / 5).astype(float)
        region = torch.tensor(region)

        k = 2 * torch.pi / self.lam
        Ny, Nx = np.shape(self.intensity_target)

        dx = self.lx / Nx
        dy = self.ly / Ny

        # convert array to work with torch
        amplitude_source = torch.as_tensor(np.sqrt(self.intensity_source / np.sum(self.intensity_source)),
                                           dtype=torch.float64)  # amplitude of the incident wave

        intensity_target = torch.tensor(self.intensity_target, requires_grad=True, dtype=torch.float64)
        phi_n = torch.tensor(self.phi_initial, requires_grad=True, dtype=torch.float64)

        # create the meshgrid
        x = torch.linspace(start=-self.lx / 2, end=self.lx / 2, steps=Nx, requires_grad=True)
        y = torch.linspace(start=-self.ly / 2, end=self.ly / 2, steps=Ny, requires_grad=True)
        xv, yv = torch.meshgrid(x, y, indexing='xy')

        # spatial frequencies
        kx = torch.as_tensor(sp.fft.fftfreq(Nx, dx) * (2 * np.pi), dtype=torch.float64)
        ky = torch.as_tensor(sp.fft.fftfreq(Ny, dy) * (2 * np.pi), dtype=torch.float64)
        kxv, kyv = torch.meshgrid(kx, ky, indexing='xy')

        impulse_response_straight = -torch.exp(torch.tensor(1j * k * self.z)) * torch.exp(
            -1j * (self.z / (2 * k)) * (torch.pow(kxv, 2) + torch.pow(kyv, 2)))

        focus = self.z
        # transmission function
        lens = torch.pow(xv, 2) + torch.pow(yv, 2) <= torch.pow(torch.tensor(50 / 2), 2)
        lens = lens.type(torch.float64)

        phi_lens = (-k / (2 * focus) * (torch.pow(xv, 2) + torch.pow(yv, 2))) % (2 * torch.pi)

        # phase function for thin lens with focal distance named focus
        phase = torch.exp(1j * phi_lens) * lens

        target_norm = torch.sum(intensity_target)
        target_normalised = intensity_target / target_norm

        def target_function(phi):
            field_n = torch.fft.ifft2((torch.fft.fft2(
                torch.fft.ifft2(torch.fft.fft2(amplitude_source * torch.exp(1j * phi)) * impulse_response_straight) * phase) *
                                   impulse_response_straight))
            intensity = torch.pow(torch.abs(field_n), 2)
            intensity_norm = torch.sum(intensity)
            intensity_normalised = intensity / intensity_norm
            return torch.sum(torch.pow((target_normalised - intensity * region), 2))

        eps_n = target_function(phi_n)

        count = 0
        dc_dphi_n = torch.autograd.grad(eps_n, phi_n)[0]
        s_n = -eps_n / torch.pow(torch.norm(dc_dphi_n, p=2), 2)

        while True:
            phi_new = phi_n + s_n * dc_dphi_n
            count += 1
            # print('C='+str(eps_n))
            if torch.abs(eps_n) <= self.accuracy:
                break
            else:

                phi_n = torch.clone(phi_new)
                eps_n = target_function(phi_new)
                dc_dphi_n = torch.autograd.grad(eps_n, phi_new)[0]
                s_n = -eps_n / torch.pow(torch.norm(dc_dphi_n, p=2), 2)
                print('iteration:'+str(count))

        phi_new_numpy = phi_new.detach().numpy() % (2 * np.pi)

        step = 2*np.pi/self.levels
        phi_binary = (phi_new_numpy//step)
        phi = phi_binary*step

        return phi, phi_binary


class GerchbergSaxtonAlgorithm:
    def __init__(self, intensity_target, intensity_source, z, lx, ly, lam, accuracy, phi_initial, levels=256):
        """
        :param intensity_target: target intensity on the screen
        :param intensity_source: intensity in the plane of the DOE
        :param z: beam spreading distance
        :param lx: x-axis length of the screen
        :param ly: y-axis length of the screen
        :param lam: wavelength
        :param accuracy: optimization termination criterion
        :param phi_initial: initial approximation of the DOE phase function
        :param levels: number of quantization levels
        """
        self.intensity_target: np.ndarray = intensity_target
        self.intensity_source: np.ndarray = intensity_source
        self.z: float = z
        self.lx: float = lx
        self.ly: float = ly
        self.lam: float = lam
        self.accuracy: float = accuracy
        self.phi_initial: np.ndarray = phi_initial
        self.levels: int = levels

    def phase_retrieval_1d(self):
        """
        optimisation algorithm for 1d phase grating

        :return: phase function and binary mask for DOE
        """

        Ny, Nx = np.shape(self.intensity_target)

        x_nf: np.ndarray = np.linspace(-self.lx / 2, self.lx / 2, Nx)
        y_nf: np.ndarray = np.linspace(-self.ly / 2, self.ly / 2, Ny)
        xv_nf, yv_nf = np.meshgrid(x_nf, y_nf)
        y = y_nf[int(Ny/2)]

        dx = self.lx / Nx
        dy = self.ly / Ny

        intensity_target_1d = self.intensity_target[int(Ny/2)]
        intensity_source_1d = self.intensity_source[int(Ny/2)]

        # spatial frequencies
        k = 2 * np.pi / self.lam
        kx_nf = sp.fft.fftfreq(Nx, dx) * (2 * np.pi)
        ky_nf = sp.fft.fftfreq(Ny, dy) * (2 * np.pi)
        kxv_nf, kyv_nf = np.meshgrid(kx_nf, ky_nf)
        ky = ky_nf[int(Ny/2)]

        amplitude_source = np.sqrt(intensity_source_1d)
        amplitude_target = np.sqrt(intensity_target_1d)

        phi_new = self.phi_initial[int(Ny/2)]

        # impulse response functions for forward and backward propagation
        impulse_response_straight = -np.exp(1j * k * self.z) * np.exp(
            -1j * (self.z / (2 * k)) * (np.power(kx_nf, 2) + np.power(ky, 2)))
        impulse_response_inverse = np.exp(-1j * k * self.z) * np.exp(
            1j * (self.z / (2 * k)) * (np.power(kx_nf, 2) + np.power(ky, 2)))

        # field in DOE plane
        E1 = amplitude_source * np.exp(1j * phi_new)

        count = 0
        errp = 1
        while True:
            # propagation from DOE to lens
            E1_fft = sp.fft.fft(E1, workers=8)
            E2_fft = E1_fft * impulse_response_straight
            E2 = sp.fft.ifft(E2_fft, workers=8)

            focus = self.z
            # field after lens
            lens = np.power(x_nf, 2) + np.power(y, 2) <= np.power(50 / 2, 2)
            lens = lens.astype(float)
            phi = (-k / (2 * focus) * (np.power(x_nf, 2) + np.power(y, 2)))*lens % (2 * np.pi)
            phase = np.exp(1j * phi)
            E3 = E2*phase

            # propagation from lens to screen
            E3_fft = sp.fft.fft(E3)
            E4_fft = E3_fft * impulse_response_straight
            E4 = sp.fft.ifft(E4_fft)

            module = np.absolute(E4)

            # new field transformation in the screen plane
            E4_new = amplitude_target*E4/module

            # backward propagation
            E4_new_fft = sp.fft.fft(E4_new)
            E3_new_fft = impulse_response_inverse * E4_new_fft
            E3_new = sp.fft.ifft(E3_new_fft)

            phi_back = (k / (2 * focus) * (np.power(x_nf, 2) + np.power(y, 2)))*lens % (2 * np.pi)
            phase_back = np.exp(1j * phi_back)

            E2_new = E3_new*phase_back

            E2_new_fft = sp.fft.fft(E2_new)
            E1_new_fft = E2_new_fft * impulse_response_inverse
            E1_new = sp.fft.ifft(E1_new_fft)

            intensity_target_new = np.power(np.abs(E4), 2)

            delta = func.relative_error(intensity=intensity_target_new, target=intensity_target_1d)
            if np.absolute(delta-errp) <= self.accuracy:
                phi_new = (np.angle(E1) + 2 * np.pi) % (2 * np.pi)
                print("Relative error(GS): "+str(delta))
                break
            else:
                print("Relative error(GS): "+str(delta))
                E1 = amplitude_source * np.exp(1j*((np.angle(E1_new) + 2*np.pi) % (2*np.pi)))
                count += 1
                errp = delta

        # phase quantisation
        phi_new = phi_new % (2*np.pi)

        step = 2 * np.pi / self.levels
        phi_1d = (phi_new // step) * step
        phi_binary_1d = (phi_new // step)

        phi = np.array([])
        phi_binary = np.array([])

        for i in range(Ny):
            phi = np.append(phi, phi_1d)
            phi_binary = np.append(phi_binary, phi_binary_1d)

        phi = np.reshape(phi, (Ny, Nx))
        phi_binary = np.reshape(phi_binary, (Ny, Nx))
        return phi, phi_binary

    def phase_retrieval_2d(self):
        """
        optimisation algorithm for 2d phase grating

        :return: phase function and binary mask for DOE
        """

        Ny, Nx = np.shape(self.intensity_target)
        dx = self.lx / Nx
        dy = self.ly / Ny

        # spatial frequencies
        k = 2 * np.pi / self.lam
        kx_nf = sp.fft.fftfreq(Nx, dx) * (2 * np.pi)
        ky_nf = sp.fft.fftfreq(Ny, dy) * (2 * np.pi)
        kxv_nf, kyv_nf = np.meshgrid(kx_nf, ky_nf)

        amplitude_source = np.sqrt(self.intensity_source)
        amplitude_target = np.sqrt(self.intensity_target)

        phi_new = self.phi_initial

        # impulse response functions for forward and backward propagation
        impulse_response_straight = -np.exp(1j * k * self.z) * np.exp(
            -1j * (self.z / (2 * k)) * (np.power(kxv_nf, 2) + np.power(kyv_nf, 2)))
        impulse_response_inverse = np.exp(-1j * k * self.z) * np.exp(
            1j * (self.z / (2 * k)) * (np.power(kxv_nf, 2) + np.power(kyv_nf, 2)))

        # field in DOE plane
        E1 = amplitude_source * np.exp(1j * phi_new)

        count = 0
        errp = 1
        while True:

            # propagation from DOE to lens
            E1_fft = sp.fft.fft2(E1, workers=8)
            E2_fft = E1_fft*impulse_response_straight
            E2 = sp.fft.ifft2(E2_fft, workers=8)

            focus = self.z
            doe_lens_forward = fo.ThinLens(50, E2, focus, self.lam, self.lx, self.ly)

            # field after the lens
            E3 = doe_lens_forward.output_field_forward()[0]

            # propagation from lens to screen
            E3_fft = sp.fft.fft2(E3)
            E4_fft = E3_fft * impulse_response_straight
            E4 = sp.fft.ifft2(E4_fft)

            module = np.absolute(E4)

            # new field transformation in the screen plane
            E4_new = amplitude_target * E4 / module

            # backward propagation
            E4_new_fft = sp.fft.fft2(E4_new)
            E3_new_fft = impulse_response_inverse * E4_new_fft
            E3_new = sp.fft.ifft2(E3_new_fft)

            doe_lens_backward = fo.ThinLens(50, E3_new, focus, self.lam, self.lx, self.ly)
            E2_new = doe_lens_backward.output_field_backward()[0]

            E2_new_fft = sp.fft.fft2(E2_new)
            E1_new_fft = E2_new_fft * impulse_response_inverse
            E1_new = sp.fft.ifft2(E1_new_fft)

            intensity_target_new = np.power(np.abs(E4), 2)

            delta = func.relative_error(intensity=intensity_target_new, target=self.intensity_target)
            if np.absolute(delta - errp) <= self.accuracy:
                phi_new = (np.angle(E1) + 2 * np.pi) % (2 * np.pi)
                print("Relative error(GS): " + str(delta))
                break
            else:
                print("Relative error(GS): " + str(delta))
                E1 = amplitude_source * np.exp(1j * ((np.angle(E1_new) + 2 * np.pi) % (2 * np.pi)))
                count += 1
                errp = delta

        # phase quantisation
        phi_new = phi_new % (2 * np.pi)

        step = 2 * np.pi / self.levels
        phi_binary = phi_new // step
        phi = phi_binary*step
        return phi, phi_binary



