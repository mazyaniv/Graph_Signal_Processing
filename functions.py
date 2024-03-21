import numpy as np
import math
import scipy.signal as ss

# def get_key_by_value(dictionary, target_value):
#     for key, value in dictionary.items():
#         if value == target_value:
#             return key

# def angles_generate(pram):
#     while True:
#         range_array = np.arange(pram.teta_range[0], pram.teta_range[1], pram.Res)
#         teta = np.random.choice(range_array[1:-1], size=pram.D, replace=False)
#         if abs(teta[0] - teta[1]) > pram.delta:
#             break
#     return np.sort(teta)[::-1]

def observ(SNR, K, A): #k=snapshots
    N = A.shape[0]
    D = A.shape[1]
    real_s = np.random.normal(0, 1 / math.sqrt(2), (D, K))
    im_s = np.random.normal(0, 1 / math.sqrt(2), (D, K))
    s = real_s + 1j * im_s
    s_samp = s.reshape(D, K)
    real_n = np.random.normal(0, (10 ** (-SNR / 20)) / math.sqrt(2), (N, K))
    im_n = np.random.normal(0, (10 ** (-SNR / 20)) / math.sqrt(2), (N, K))
    n = real_n + 1j * im_n
    n_samp = n.reshape(N, K)
    x_a_samp = (A@s_samp) + n_samp
    return x_a_samp.T.flatten()

def GFT(S, x):
    eigenvalues, eigenvectors = np.linalg.eig(S)
    sorted_eigenvectors = eigenvectors[:, np.argsort(np.abs(eigenvalues))] #sorted in case of negative eigenvalues
    x_tag = sorted_eigenvectors.T@x

    return x_tag


