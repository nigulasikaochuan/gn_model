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

    space = 35e9
    power_lin = (10**(2/10))/1000

    snrs = []
    bers=[]
    for span_number in range(10,70):
        print(span_number)
        signals = [
            Signal(signal=power_lin, nli=0, ase=0, mf='dp-qpsk', carri=193.1e12 + i * space, baudrate=32e9,
                   number=i)
            for i in range(15)]
        spans = [Span(alpha=0.2,D=16.7,gamma=1.3,length=120,lam=0) for i in range(span_number)]

        fiber = Fiber(spans=spans,edfa=Edfa(gain=spans[0].length*spans[0].alpha,nf=5))
        fiber.prop(signals)
        center = signals[7]
        snr = center.signal/(center.ase+center.nli)
        snrs.append(snr)
        bers.append(calc_ber(snr))

    snr_db = 10*np.log10(np.array(snrs))
    # plt.semilogy(snr_db,bers)
    # plt.show()
    # plt.grid()
    plt.plot(bers)
    plt.show()