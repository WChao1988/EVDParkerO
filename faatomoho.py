from scipy.fftpack import fft2, ifft2
import numpy as np
from abc import ABCMeta, abstractmethod
from math import factorial
from scipy.signal.windows import tukey
import os


file_path = os.path.dirname(os.path.abspath(__file__))


class Faa2MultiInterface(metaclass=ABCMeta):
    def __init__(self, faa_matrix, delta_bnds, delta_rhos, mus, reference_depths, longrkm, longckm):
        self.faa_matrix = faa_matrix
        self.delta_bnds = delta_bnds
        self.delta_rhos = delta_rhos
        self.mus = mus
        self.reference_depths = reference_depths
        self.longrkm, self.longckm = longrkm, longckm

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


class EVDMOHO(Faa2MultiInterface):
    def __init__(self, faa_matrix, delta_bnds, delta_rhos, mus, reference_depths, longrkm, longckm, wh, alpha):
        super(EVDMOHO, self).__init__(faa_matrix, delta_bnds, delta_rhos, mus, reference_depths, longrkm, longckm)
        self.G = 6.67
        self.wh = wh
        self.alpha = alpha
        self.frequency = self.__frequency__()
        self.filter = self.__filter__()
        self.delta_bnds = delta_bnds

    def __frequency__(self):
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
        bnd_twkey = self.twkey(self.delta_bnds[k])
        bnd_twkey_n = bnd_twkey ** n
        bnd_fourier = fft2(bnd_twkey_n)
        return bnd_fourier

    def __faa_fft__(self):
        faa_fft = fft2(self.twkey(self.faa_matrix))
        return faa_fft

    def once_downward(self, t):
        fft = self.__faa_fft__() / 2 / np.pi / self.G
        summary0 = 0 + 0j
        for n in range(2, t + 1):
            summary0 = summary0 + (self.frequency - self.mus[0])**(n-1) / factorial(n) * self.bnd_k_n(0, n)
        summary1 = 0 + 0j
        for n in range(1, t + 1):
            summary1 = summary1 + (self.frequency - self.mus[1])**(n-1) / factorial(n) * self.bnd_k_n(1, n)
        residual_fft = (fft * np.exp(self.frequency*self.reference_depths[0]) / self.delta_rhos[0] -
                        self.delta_rhos[0]*summary0 * np.exp(self.frequency*self.reference_depths[0]-self.frequency*self.reference_depths[0]) / self.delta_rhos[0] -
                        self.delta_rhos[1]*summary1 * np.exp(self.frequency*self.reference_depths[0]-self.frequency*self.reference_depths[1]) / self.delta_rhos[0])
        residual_fft_0 = residual_fft * self.filter
        residual_fft_0[0, 0] = 0
        self.delta_bnds[0] = ifft2(residual_fft_0).real

    def downward(self, t=8):
        for _ in range(1, t+1):
            self.once_downward(_)
        return self.delta_bnds[0] - self.reference_depths[0]