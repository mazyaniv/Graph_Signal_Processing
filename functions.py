import math
import scipy.signal as ss
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from classes import prameters_class, Matrix_class

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
def generate_qpsk_symbols(K,D):
    random_symbols = np.random.randint(0, 4, (K,D))
    qpsk_constellation = np.array([1+1j, -1+1j, -1-1j, 1-1j]) / np.sqrt(2)
    qpsk_symbols = qpsk_constellation[random_symbols]
    return qpsk_symbols.T
def observ(SNR, K, A): #k=snapshots
    N = A.shape[0]
    D = A.shape[1]
    # sample_rate = 1e6
    # t = np.arange(K) / sample_rate  # time vector
    # f_tone = np.array([0.02e3,0.08e3])
    # s = np.exp(2j * np.pi * f_tone.reshape(2,1) * t.reshape(1,800))
    # real_s = np.random.normal(1, 1 / math.sqrt(2), (D, K))
    # im_s = np.random.normal(1, 1 / math.sqrt(2), (D, K))
    # s = real_s + 1j * im_s
    s = generate_qpsk_symbols(K,D)
    s_samp = s.reshape(D, K)

    real_n = np.random.normal(0, (10 ** (-SNR / 20)) / math.sqrt(2), (N, K))
    im_n = np.random.normal(0, (10 ** (-SNR / 20)) / math.sqrt(2), (N, K))
    n = real_n + 1j * im_n
    n_samp = n.reshape(N, K)
    x_a_samp = (A@s_samp) + n_samp
    return x_a_samp

def GFT(S, x):
    eigenvalues, eigenvectors = np.linalg.eig(S)
    sorted_eigenvectors = eigenvectors[:, np.argsort(np.real(eigenvalues))] #sorted in case of negative eigenvalues
    # print(np.imag(eigenvalues))
    # print(np.argsort(np.abs(eigenvalues)))
    x_tag = sorted_eigenvectors.T.conjugate()@x
    vector = np.zeros(len(x_tag), dtype=int)
    vector[:] = 1
    B = np.diag(vector)
    x_tag = B@x_tag
    return x_tag

def quantize(A, P, thresh_real=0, thresh_im=0):
    mask = np.zeros(np.shape(A), dtype=complex)
    mask[:P, :] = (1 / math.sqrt(2)) * (np.sign(A[:P, :].real - (thresh_real)) + (1j * (np.sign(A[:P, :].imag - ((thresh_im))))))
    mask[P:, :] = A[P:, :]
    return mask
def G_DOA(pram,teta_range,S,q):
    theta_range = np.arange(teta_range[0], teta_range[1], pram.Res)
    A_s_mat = np.zeros((pram.M, pram.M,len(theta_range)), dtype=complex)
    for theta in theta_range:
        my_parameters2 = prameters_class(pram.M, pram.SNR, pram.K, [theta])
        steering2 = Matrix_class(my_parameters2).steering()
        A_s = steering2 @ steering2.T.conjugate() - (1 / (pram.M - 1)) * np.eye(pram.M)  # Matrix_class(my_parameters2).Adjacency()
        A_s_mat[:,:,theta] = A_s
    spectrum_vec = np.zeros((pram.K,len(theta_range)))
    obs_a = observ(pram.SNR, pram.K, S)
    x_vec = quantize(obs_a, q)
    for i in range(pram.K):
        spectrum = np.zeros(len(theta_range))
        for idx, theta in enumerate(theta_range):
            x_tag_s = GFT(A_s_mat[:,:,theta], x_vec[:,i])
            sorted_indices = np.argsort(x_tag_s)[::-1]
            # spectrum[idx]= 1/LA.norm(np.abs(x_tag_s[:-pram.D])/np.max(np.abs(x_tag_s)))
            spectrum[idx]= 1/LA.norm(np.abs(np.delete(x_tag_s, sorted_indices[:1]))/np.max(np.abs(x_tag_s)))
        spectrum_vec[i,:]=spectrum
    spectrum = np.average(spectrum_vec,0)
    plt.plot(theta_range, spectrum,marker="*")
    plt.show()
    peaks, _ = ss.find_peaks(spectrum,distance=pram.delta/pram.Res)
    peaks = list(peaks)
    peaks.sort(key=lambda x: spectrum[x])
    pred = np.array(peaks[-pram.D:])
    pred = np.sort(pred)[::-1]
    return pred*pram.Res+teta_range[0]


