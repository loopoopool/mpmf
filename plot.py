import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('iei.dat').transpose()

k1 = 3
k2 = 3

tmp = np.where( data[0] == k1, data, np.zeros(data.shape))
tmp = np.where( tmp[1] == k2, tmp, np.zeros(data.shape) )

zero_columns_mask = np.all(tmp == 0, axis=0)
tmp = tmp[:, ~zero_columns_mask]
tmp = tmp[2:]

mat = np.zeros((2*k1+1,2*k2+1))
for i, x in enumerate(tmp[2]):
        mat[int(tmp[0, i])+k1, int(tmp[1, i])+k2] = x

plt.imshow(mat, cmap='seismic', interpolation='nearest', vmin=-3, vmax=3)
plt.colorbar()
plt.show()
