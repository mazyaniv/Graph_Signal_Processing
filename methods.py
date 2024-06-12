import numpy as np
from matplotlib import pyplot as plt
import math
import scipy.signal as ss
def root_music(pram,R): # Root-MUSIC
    my_vec_coff = np.zeros((pram.M, 2 * pram.M - 1), dtype=complex)
    eigvals, eigvecs = np.linalg.eig(R)
    sorted_indices = np.argsort(eigvals.real)[::-1]  # Sort eigenvalues in descending order
    eigvecs_sorted = eigvecs[:, sorted_indices]
    En = eigvecs_sorted[:, pram.D:]
    matrix = En@En.conj().T
    for i in range(np.shape(matrix)[1]):
        vector = np.concatenate((matrix[:, i][::-1].reshape(1, -1), np.zeros((1, pram.M - 1))), axis=1)
        my_vec_coff[i, :] = np.roll(vector, i)

    cofficients = np.sum(my_vec_coff, 0)
    roots = np.poly1d(cofficients[::-1]).r
    sorted_roots = sorted(roots, key=lambda r: abs(abs(r) - 1))
    closest_roots = sorted_roots[0], sorted_roots[2]  # two closest roots
    pred = -np.degrees(np.arcsin(np.angle(closest_roots) / math.pi))[::-1]  # TODO why "-"?
    pred = np.sort(pred)[::-1]
    return pred
def esprit(pram,R):
    eigvals, eigvecs = np.linalg.eig(R)
    sorted_indices = np.argsort(eigvals.real)[::-1]  # Sort eigenvalues in descending order
    eigvecs_sorted = eigvecs[:, sorted_indices]
    Es = eigvecs_sorted[:, :pram.D]
    S1 = Es[1:,:]
    S2 = Es[:-1,:]
    P = np.linalg.inv(S1.conj().transpose()@S1)@S1.conj().transpose()@S2 #LS
    eigvals, eigvecs = np.linalg.eig(P)
    pred = np.degrees(np.arcsin(np.angle(eigvals) / math.pi))
    pred = np.sort(pred)[::-1]
    return pred
def music(pram,R):
    eigvals, eigvecs = np.linalg.eig(R)
    sorted_indices = np.argsort(eigvals.real)[::-1]  # Sort eigenvalues in descending order
    eigvecs_sorted = eigvecs[:, sorted_indices]
    En = eigvecs_sorted[:, pram.D:]

    theta_range = np.radians(np.arange(pram.teta_range[0], pram.teta_range[1], pram.Res))  # Convert angles to radians
    music_spectrum = np.zeros(len(theta_range))
    for idx, theta in enumerate(theta_range):
        steering_vector = np.exp(-1j * np.pi * np.arange(pram.M) * np.sin(theta))
        music_spectrum[idx] = 1 / np.linalg.norm(En.conj().T @ steering_vector)
    # plt.plot(np.degrees(theta_range), music_spectrum)
    # plt.show()

    peaks, _ = ss.find_peaks(music_spectrum)
    peaks = list(peaks)
    peaks.sort(key=lambda x: music_spectrum[x])
    pred = np.array(peaks[-pram.D:])
    pred = np.sort(pred)[::-1]
    return pred*pram.Res+pram.teta_range[0]