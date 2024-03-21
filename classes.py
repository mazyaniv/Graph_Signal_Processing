import numpy as np

class prameters_class():
    def __init__(self, M,N_q,SNR,snapshot,D,teta_range,monte,delta,Res,dictio):
        self.M = M
        self.N_q = N_q
        self.D = D
        self.teta_range = teta_range
        self.SNR = SNR
        self.snapshot = snapshot
        self.monte = monte
        self.delta = delta
        self.Res = Res
        self.dictio = dictio
class Matrix_class():
    def __init__(self, N, K, omega_0, tau):
        self.N = N
        self.K = K
        self.omega_0 = omega_0
        self.tau = tau
        self.matrix = np.zeros((N, N), dtype=complex)
        self.matrix2 = np.zeros((K, K), dtype=complex)
    def generate_matrix(self):
        for i in range(self.N):
            if i == 0:
                self.matrix[i, 1] = np.exp(1j * self.omega_0 * self.tau)
                self.matrix[i, -1] = np.exp(1j * self.omega_0 * self.tau * (self.N - 1))
            elif i == self.N - 1:
                self.matrix[i, -2] = np.exp(-1j * self.omega_0 * self.tau)
                self.matrix[i, 0] = np.exp(-1j * self.omega_0 * self.tau * (self.N - 1))
            else:
                self.matrix[i, i - 1] = np.exp(-1j * self.omega_0 * self.tau)
                self.matrix[i, i + 1] = np.exp(1j * self.omega_0 * self.tau)
        A_time = 0.5 * self.matrix

        for i in range(self.K):
            if i == 0:
                self.matrix2[i, 1] = np.exp(-1j * self.omega_0)
                self.matrix2[i, -1] = np.exp(-1j * self.omega_0 * (self.K - 1))
            elif i == self.K - 1:
                self.matrix2[i, -2] = np.exp(1j * self.omega_0)
                self.matrix2[i, 0] = np.exp(1j * self.omega_0 * (self.K - 1))
            else:
                self.matrix2[i, i - 1] = np.exp(1j * self.omega_0)
                self.matrix2[i, i + 1] = np.exp(-1j * self.omega_0)
        A_space = 0.5 * self.matrix2
        return np.kron(A_time, A_space)

if __name__ == "__main__":
    print("Not main file")