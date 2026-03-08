from scipy.fftpack import fft2, ifft2
import numpy as np
from abc import ABCMeta, abstractmethod
from math import factorial
from scipy.signal.windows import tukey
import os


file_path = os.path.dirname(os.path.abspath(__file__))


class FaaToMultiInterface(metaclass=ABCMeta):
    def __init__(self, faa_matrix, sediment_matrix, delta_rhos, delta_rho_initial, mus, reference_depths, longrkm, longckm, longwave_topo):
        """
        initialize the inputs and parameters
        :param faa: A matrix for gravity anomalies.
        :param delta_bnds: A dictionary that stores the undulations of the density interfaces
                           except for the density interface to be inverted.
        :param delta_rhos: A dictionary that stores the density contrasts for the density interfaces.
        :param reference_depths: A dictionary that stores the reference depths of the density interfaces.
        :param longrkm: float for row length
        :param longckm: float for col length
        """
        self.faa_matrix = faa_matrix
        self.sediment_matrix = sediment_matrix
        self.delta_rhos = delta_rhos
        self.reference_depths = reference_depths
        self.longrkm, self.longckm = longrkm, longckm
        self.mus = mus
        self.delta_rho_initial = delta_rho_initial
        self.longwave_topo = longwave_topo

    @abstractmethod
    def downward(self, t):
        """
        calculate the density interface of interest with downward iteration steps
        :param t: iteration for downwards
        :param criteria: criteria for downwards iteration
        :return: matrix for undulation of density interface of interest
        """
        pass

    @classmethod
    def twkey(cls, matrix, edge=0.02):
        """
        tukey the border values for smoothness with outside values 0
        :param matrix: matrix for tukey
        :param edge: float for definition of border
        :return: tukey matrix
        """
        nrow, ncol = matrix.shape
        tky = np.array([row_tky * col_tky for row_tky in tukey(nrow, edge)
                        for col_tky in tukey(ncol, edge)]).reshape(matrix.shape)
        return tky * matrix


class EVDBATHYMETRY(FaaToMultiInterface):
    def __init__(self, faa_matrix, sediment_matrix, delta_rhos, delta_rho_initial, mus, reference_depths, longrkm, longckm, wh, alpha, longwave_topo):
        """
        :param faa_matrix: a matrix for faa
        :param sediment_matrix:  a matrix for sediment
        :param delta_rhos: a dict for density contrasts
        :param delta_rho_initial: a float for the density contrast with the initialization
        :param mus: a float for factor of decrease
        :param reference_depths: a dict for reference depths
        :param longrkm:
        :param longckm:
        :param wh: a float for cutoff of high frequency
        :param alpha: an int for the compress of the high frequency
        """
        super(EVDBATHYMETRY, self).__init__(faa_matrix, sediment_matrix, delta_rhos, delta_rho_initial, mus, reference_depths, longrkm, longckm, longwave_topo)
        self.G = 6.67
        self.wh = wh
        self.alpha = alpha
        self.frequency = self.__frequency__()
        self.filter = self.__filter__()
        self.delta_bnds = {}
        # temp
        self.temp = {}

    def __frequency__(self):
        """
        inner function for calculating the frequency
        :return: frequency matrix
        """
        nrow, ncol = self.faa_matrix.shape
        frequency = np.zeros((nrow, ncol))
        for i in range(nrow):
            for j in range(ncol):
                ii = i if i <= nrow / 2 else i - nrow
                jj = j if j <= ncol / 2 else j - ncol
                frequency[i, j] = 2 * np.pi * np.sqrt((ii / self.longrkm) ** 2 + (jj / self.longckm) ** 2)
        return frequency

    def __filter__(self):
        """
        inner function for calculating the filter that lowpass the frequency value
        :return: a filter matrix
        """
        nrow, ncol = self.faa_matrix.shape
        filter = np.ones(self.faa_matrix.shape)
        for i in range(nrow):
            for j in range(ncol):
                if self.frequency[i, j] > self.wh:
                    ratio = self.frequency[i, j] / self.wh
                    filter[i, j] = ratio ** (1 - self.alpha) - (1 - self.alpha) * np.log(ratio) * ratio ** (
                                1 - self.alpha)
        return filter

    def bnd_k_n(self, k, n):
        name = '%d-%d' % (k, n)
        if self.temp.get(name) is None:
            bnd_twkey = self.twkey(self.delta_bnds[k])
            bnd_twkey_n = bnd_twkey ** n
            bnd_fourier = fft2(bnd_twkey_n)
            self.temp[name] = bnd_fourier
        else:
            bnd_fourier = self.temp.get(name)
        return bnd_fourier

    def __faa_fft__(self):
        name = 'faa-ft'
        if self.temp.get(name) is None:
            faa_ft = fft2(self.twkey(self.faa_matrix))
            self.temp[name] = faa_ft
        else:
            faa_ft = self.temp.get(name)
        return faa_ft

    def __initialize__interface__(self):
        """
        initialize the target density interface with only the gravity anomalies
        :return: matrix for interface
        """
        qg_target = 2 * np.pi * self.G * self.delta_rho_initial * np.exp(-self.frequency *
                                                                         self.reference_depths[1])
        # faa fft
        faa_fft = self.__faa_fft__()
        # interface_fft
        interface_fft = faa_fft / qg_target * self.filter
        # set the constant term been zero
        interface_fft[0, 0] = 0 + 0j
        # inverse interface
        interface1 = ifft2(interface_fft).real - self.reference_depths.get(1)
        # form the sedimentary layer
        interface0 = interface1 - self.sediment_matrix
        self.delta_bnds[0] = interface0 + self.reference_depths.get(0)
        self.delta_bnds[1] = interface1 + self.reference_depths.get(1)
        return self.delta_bnds[1]

    def once_downward(self, t):
        fft = self.__faa_fft__() / 2 / np.pi / self.G / self.delta_rhos[1] / \
              np.exp(-self.frequency * self.reference_depths[1])

        factor_0 = self.delta_rhos[0] / self.delta_rhos[1] * np.exp(self.frequency
                                                                    * (self.reference_depths[1] -
                                                                       self.reference_depths[0]))
        summary0 = 0 + 0j
        for n in range(1, t + 1):
            summary0 = summary0 + (self.frequency - self.mus[0])**(n-1) / factorial(n) * self.bnd_k_n(0, n)
        summary1 = 0 + 0j
        for n in range(2, t + 1):
            summary1 = summary1 + (self.frequency - self.mus[1])**(n-1) / factorial(n) * self.bnd_k_n(1, n)
        residual_fft = (fft - summary0 * factor_0 - summary1) * self.filter
        residual_fft[0, 0] = 0
        self.delta_bnds[1] = ifft2(residual_fft).real
        self.delta_bnds[0] = self.delta_bnds[1] - self.reference_depths[1] - \
                                                                         self.sediment_matrix + self.reference_depths[0]

    def downward(self, t=8):
        self.__initialize__interface__()
        for _ in range(2, t+1):
            self.once_downward(_)
        return self.delta_bnds[1] - self.reference_depths[1]+ self.longwave_topo


