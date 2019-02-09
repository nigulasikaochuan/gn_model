from gn_model import Span
from gn_model import Edfa
from gn_model import Signal
from gn_model import Fiber


import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    power = (10**(2/10))/1000
    space = 33.6e9
    #sig = [Signal(signal=power,nli=0,ase=0,space=space,carri=193.1e12+i*space,baudrate=32e9,number=i) for i in range(9)]


    eta = []

    edfa = Edfa(gain=20, nf=5)
    for span_index in range(50):
        sig = [Signal(signal=power, nli=0, ase=0, space=space, carri=193.1e12 + i * space, baudrate=32e9, number=i) for
               i in range(9)]

        span = [Span(alpha=0.2,D=16.7,gamma=1.3,length=100,lam=0) for i in range(span_index+1)]

        fiber = Fiber(spans=span,edfa=edfa)

        fiber.prop(sig)
        eta.append(10*np.log10(sig[7].nli/sig[7].signal**3))

    plt.semilogx(np.arange(50)+1,eta)
    plt.xlim(1,50)
    plt.ylim(20,50)
    plt.grid()
    plt.show()
