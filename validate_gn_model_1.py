from gn_model import Fiber
from gn_model import Signal
from gn_model import Span
from gn_model import Edfa
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc


def calc_ber(snr):
    return 0.5 * erfc(np.sqrt(0.5 * snr))


if __name__ == '__main__':

    for spa in [33.6,]
    spa =50* 1e9
    sr = []

    for power in range(-10, 5):
        num = 0
        power_lin = (10 ** (power / 10)) * 0.001
        signals = [
            Signal(signal=power_lin, nli=0, ase=0, mf='dp-qpsk', carri=193.1e12 + i * spa, baudrate=32e9, number=i)
            for i in range(9)]
        span = Span(alpha=0.2, D=16.7, gamma=1.5, length=120, lam=1550e-9)
        edfa = Edfa(gain=120 * 0.2, nf=5)
        center = signals[3]
        while True:
            span.prop(center,signals)
            edfa.prop(center)
            snr = center.signal/(center.nli+center.ase)

            if calc_ber(snr)>1e-3:
                break
            else:
                num+=1

        print(num)
