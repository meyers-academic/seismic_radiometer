from gwpy.spectrum import Spectrum
import numpy as np


class Spec(Spectrum):
    def whiten(self, width=1):
        """
        whitens spectrum by getting a smoothed version of the
        absolute value of the spectrum and dividing out by that smoothed
        version.

        Parameters
        ----------
        width : `int`, kwarg, optional, default=1
            width over which to smooth (in Hz)

        Returns
        -------
        out : `Spec`, whitened spectrum
        """

        S = self
        env = S.smooth(width=width)
        return S / env

    def smooth(self, width=1):
        """
        Smooths spectrum by convolving data with ones of
        proper length. Averages over width

        Parameters
        ----------
        width : `int`, optional, default=1,
            Number of seconds or Hz to use in convolution.

        Returns
        -------
        smoothed : `Trace`
            Smoothed time series trace.
        """

        S = np.abs(self)

        # turn width into # of samples
        width = width * (1 / S.df.value)

        # get window
        window = np.ones((2 * width))

        # do convolution
        S = np.convolve(S, window / (2 * width), 'same')
        S = Spec(S)
        S.__dict__ = self.copy_metadata()
        return S
