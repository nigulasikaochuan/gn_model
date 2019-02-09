import matplotlib.pyplot as plt

snr_gn = []
with open('snr') as snrs:
    for snr in snrs:
        snr_gn.append(float(snr))

plt.plot(range(-5,5),snr_gn,marker='o',label='GN Model')
plt.legend()
plt.show()
plt.xlabel('$launch power [dbm]')
plt.xlim(-5,5)
plt.ylabel('snr[db]')
