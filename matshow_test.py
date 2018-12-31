import matplotlib.pyplot as plt
import numpy as np


def samplemat(dims):
    """Make a matrix with all zeros and increasing elements on the diagonal"""
    aa = np.zeros(dims)
    for i in range(min(dims)):
        aa[i, i] = i**3
    return aa


# Display matrix
print(samplemat((15,15)))
plt.matshow(samplemat((15, 15)))

plt.show()
