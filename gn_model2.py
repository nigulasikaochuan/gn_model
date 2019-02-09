import numpy as np
from scipy.constants import c
from scipy.constants import h
# c = 3e8
class Carriers(object):

    def __init__(self,nch,df,baud_rate,power,center_frequence,roll_off):
        self.nch = nch
        self.baud_rate = baud_rate
        self.df = df
        self.power = power
        self.center_frequence = center_frequence
        self.roll_off = roll_off

    @property
    def gwdm(self):
        return self.nch*self.power/self.bwdm

    @property
    def baud_rate_in_hz(self):
        return self.baud_rate

    @property
    def bwdm(self):
        return self.nch *self.bch
    @property
    def bch(self):
        return (1+self.roll_off)*self.baud_rate_in_hz
class Fiber(object):

    def __init__(self,alpha,D,gamma,length):
        self.alpha = alpha #db/km
        self.D = D
        self.gamma = gamma #  [W^-1*km^-1]
        self.length = length #km

    @property
    def alpha_lin(self):

        return  self.alpha / (2 * 10*np.log10(np.exp(1)))
    @property
    def leff(self):
        return (1-np.exp(-2*self.alpha_lin*self.length))/2/self.alpha_lin

    @property
    def leff_a(self):
        return 1/2/self.alpha_lin

    @property
    def beta2(self):
        '''

        :return: 二阶色散
        '''
        return -self.D * (1550e-9 * 1e-3) ** 2 / (2 * np.pi * c * 1e-3)  # s**2/km

    def gnli_center(self,carriers:Carriers):

        gnli = self.gamma**2*carriers.gwdm**3*self.leff**2
        gnli = gnli*np.arcsinh(0.5*np.pi**2*abs(self.beta2)*self.leff_a*carriers.bch**2*carriers.nch**(2*carriers.bch/carriers.df))
        gnli = gnli/(np.pi*abs(self.beta2)*self.leff_a)
        return 8/27*gnli

    def gnli_noise(self,carriers):
        return self.gnli_center(carriers)*carriers.bch

def calu_ase(mu,gain,nf,baud_rate):
    nf_lin = 10**(nf/10)
    gain_lin = 10**(gain/10)
    s = h*mu*(gain_lin-1)*(nf_lin*gain_lin-1)/(2*(gain_lin-1))

    return 2*s*baud_rate

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    for power in range(-6,1):

        power_lin = (10**(power/10))/1000
        span_number = 10
        carriers = Carriers(nch=11,df=50e9,baud_rate=34.5e9,power=power_lin,center_frequence=c/1550.12e-9,roll_off=0.14)


        span = Fiber(alpha=0.21,D=2.8,gamma=2,length=80)

        noise = span.gnli_noise(carriers)*span_number

        noise = noise+calu_ase(carriers.center_frequence,16,5,carriers.baud_rate)*0

        snr = power_lin/noise
        print(10*np.log10(snr))
