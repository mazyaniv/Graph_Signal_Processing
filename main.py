import numpy as np

def generate_matrix_time(N,K, omega_0, tau):
    matrix = np.zeros((N, N), dtype=complex)
    for i in range(N):
        if i == 0:
            matrix[i, 1] = np.exp(1j * omega_0 * tau)
            matrix[i, -1] = np.exp(1j * omega_0 * tau * (N - 1))
        elif i == N - 1:
            matrix[i, -2] = np.exp(-1j * omega_0 * tau)
            matrix[i, 0] = np.exp(-1j * omega_0 * tau * (N - 1))
        else:
            matrix[i, i - 1] = np.exp(-1j * omega_0 * tau)
            matrix[i, i + 1] = np.exp(1j * omega_0 * tau)
    A_time = 0.5*matrix

    matrix2 = np.zeros((K, K), dtype=complex)
    for i in range(K):
        if i == 0:
            matrix2[i, 1] = np.exp(-1j * omega_0)
            matrix2[i, -1] = np.exp(-1j * omega_0 * (K - 1))
        elif i == K - 1:
            matrix2[i, -2] = np.exp(1j * omega_0)
            matrix2[i, 0] = np.exp(1j * omega_0 * (K - 1))
        else:
            matrix2[i, i - 1] = np.exp(1j * omega_0)
            matrix2[i, i + 1] = np.exp(-1j * omega_0)
    A_space =  0.5 * matrix2
    return np.kron(A_time, A_space)


