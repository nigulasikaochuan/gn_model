import numpy as np
from numpy import pi
from scipy.constants import h, c
from numpy import arcsinh
from typing import List

import matplotlib.pyplot as plt


class Signal(object):

    '''
        Signal Class:
        property:
            signal : the launch power of signal
            nli: the non linear noise of signal
            ase: the ase noise of edfa

            carri: the carrier frequence
            baud_rate : symbol rate
            number: channel number
            mf: modulation-format
    '''

    def __init__(self, **kwargs):
        self.signal = kwargs.get('signal', 0)  # w
        self.nli = kwargs.get('nli', 0)  # w
        self.ase = kwargs.get('ase', 0)  # w
        self.carri = kwargs.get('carri', 193.1e12)  # hz
        self.baud_rate = kwargs.get('baudrate', 35e9)
        self.number = kwargs.get('number', 0)
        self.mf = kwargs.get('mf', 'dp-qpsk')
        self.roll_off = 0.02

    @property
    def bch(self):
        return self.baud_rate*(1+self.roll_off)

    @property
    def psd(self):
        '''

        :return:信号的功率谱密度， 信号功率/信号带宽
        '''
        return self.signal / self.bch  # w/hz

    @property
    def ase_psd(self):
        return self.ase/self.baud_rate

    def __str__(self):
        stra = f' signal power is {self.signal}w,nli power is {self.nli}w, ase power is {self.ase}w, ' \
               f'carrier_frequence is {self.carri} ' \
               f'baud_rate is {self.baud_rate},number is {self.number},mf is {self.mf}'
        return stra


class Span(object):

    '''
        Span class:
        property:
            length:span's length
            D: dispersion
            gamma: nonlinear
            lam: signal wavelength
            alpha:
    '''
    def __init__(self, **kwargs):
        self.length = kwargs.get('length', 80)  # km
        self.D = kwargs.get('D', 17)  # [ps/nm/km]
        self.gamma = kwargs.get('gamma', 1.2)  # [W^-1*km^-1]
        self.lam = kwargs.get('lam', 1550e-9)  # signal wavelength
        self.alpha = kwargs.get('alpha', 0.2)  # alpha db/km

    @property
    def leff(self):
        '''

        :return:有效长度
        '''
        return (1 - np.exp(-2 * self.alpha_lin * self.length)) / 2 / self.alpha_lin

    @property
    def beta2(self):
        '''

        :return: 二阶色散
        '''
        return -self.D * (1550e-9 * 1e-3) ** 2 / (2 * pi * c * 1e-3)  # s**2/km

    @property
    def alpha_lin(self):
        '''

        :return:得到线性的衰减系数
        '''
        return np.log(10 ** (self.alpha / 20))  # 1/km

    def linear_prop(self, signal: Signal):
        '''

        :param signal: Signal object
        :return: None

        线性传输，也就是功率衰减
        '''

        signal.signal = signal.signal * np.exp(-2 * self.alpha_lin * self.length)
        signal.ase = signal.ase * np.exp(-2 * self.alpha_lin * self.length)

    def prop(self, signal: Signal, signals: List[Signal]):
        '''

        :param signal:Signal object in study
        :param signals: all wdm signal
        :return: 非线性噪声

        先计算出非线性噪声，计算出的是接收机处看到的非线性噪声
        然后信号功率，ase噪声衰减
        '''
        # self.lam = c / signal.carri
        gnli = 0
        for interfer_signal in signals:
            if interfer_signal.number == signal.number:
                phi = self.calc_phi_self(interfer_signal)
                gnli += interfer_signal.psd ** 3 * phi
            else:
                phi = self.calc_phi_other(signal, interfer_signal)
                gnli += interfer_signal.psd ** 2 * signal.psd * phi
        gnli *= (16 / 27) * self.gamma * self.gamma * self.leff * self.leff

        nli_noise = gnli * signal.baud_rate
        signal.nli = nli_noise + signal.nli
        self.linear_prop(signal)


        return nli_noise

    def calc_phi_self(self, intf: Signal):

        fenzi = arcsinh(pi * pi / 2 * abs(self.beta2) * ((2 * self.alpha_lin) ** (-1)) * intf.bch ** 2)
        fenmu = 2 * pi * abs(self.beta2) * (2 * self.alpha_lin) ** (-1)
        return fenzi / fenmu

    def calc_phi_other(self, signal: Signal, interfer_signal: Signal):
        diyixiang_fenzi = arcsinh(pi ** 2 * 0.5 * self.alpha_lin ** (-1) * abs(self.beta2) * (
                interfer_signal.carri - signal.carri + interfer_signal.bch / 2) * signal.bch)
        fenmu = 4 * pi * (2 * self.alpha_lin) ** (-1) * abs(self.beta2)

        diyixiang = diyixiang_fenzi / fenmu

        dierxiang_fenzi = arcsinh((pi * pi) * (2 * self.alpha_lin) ** (-1) * abs(self.beta2) * (
                interfer_signal.carri - signal.carri - interfer_signal.bch / 2) * signal.bch)
        dierxiang = dierxiang_fenzi / fenmu

        return 2 * (diyixiang - dierxiang)


class Edfa(object):

    def __init__(self, gain, nf):
        self.gain = gain
        self.nf = nf

    @property
    def gain_lin(self):
        return 10 ** (self.gain / 10)


    def ase(self, signal):
        nf_lin = 10 ** (self.nf / 10)
        gain_lin = 10 ** (self.gain / 10)
        s = h * (c/1550e-9) * (nf_lin * gain_lin - 1) / (2 )

        # return 2*s
        return 0

    def prop(self, signal: Signal):
        '''
        因为非线性噪声是在接收机处看到的噪声，所以不考虑非线性噪声的放大和衰减
        :param signal:
        :return:
        '''
        signal.ase = signal.ase * self.gain_lin  # 之前EDFA产生的ase先被放大
        signal.ase += self.ase(signal) * signal.bch
        signal.signal = signal.signal * self.gain_lin  # 信号功率放大

        # signal.nli = signal.nli*self.gain_lin

class Fiber(object):

    def __init__(self, spans: List[Span], edfa: Edfa):
        '''

        :param spans:List of span object
        :param edfa: assume edfa is the same after every span
        '''

        self.spans = spans
        self.edfa = edfa

    def prop(self, signals: List[Signal]):
        '''

        :param signals:所有的wdm信道
        :return: None
        对于所有的wdm信号中的每一个信号，以次通过所有的span
        '''
        for signal in signals:
            for span in self.spans:
                span.prop(signal, signals)
                self.edfa.prop(signal)


if __name__ == '__main__':

    snr = []

    space = 75e9
    for power in range(0,1):

        power_lin = (10 ** (power/ 10)) * 0.001




        sigs = [
            Signal(signal=power_lin, nli=0, ase=0, carri=193.1e12 + j * space, baudrate=70e9, number=j, mf='dp-16qam')

            for j in range(9)]

        spans = [Span(length=80, D=16.7, gamma=1.3, lam=1550e-9, alpha=0.2) for k in range(15)]
        edfa = Edfa(gain=16, nf=5)
        center_channel = sigs[4]

        for span in spans:
            span.prop(center_channel,sigs)
            edfa.prop(center_channel)

        snr.append(center_channel.signal/(center_channel.nli+0))

    print(10*np.log10(np.array(snr)))
