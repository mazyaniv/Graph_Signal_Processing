import numpy as np
import math
class prameters_class():
    def __init__(self, N,SNR,k,D,teta_range,monte,delta,Res):
        self.N = N
        self.D = D
        self.teta_range = teta_range
        self.SNR = SNR
        self.k = k
        self.monte = monte
        self.delta = delta
        self.Res = Res
class Matrix_class(): #TODO D>1
    def __init__(self, N, K,teta):
        self.N = N #sensors
        self.K = K #snapshots
        self.teta = np.radians(teta)
        self.D = len(teta) #sources
        self.Stee = np.zeros((self.N, self.D), dtype=complex) #steering matrix
        self.matrix = np.zeros((N, N), dtype=complex)
        self.matrix2 = np.zeros((K, K), dtype=complex)
    def steering(self):
        A_mask = np.zeros((self.N, self.D), dtype=complex)
        for j in range(self.D):
            A_mask[:, j] = np.exp(-1j * np.pi * np.arange(self.N) * np.sin(self.teta[j]))
        self.Stee = A_mask
        return self.Stee
    def Adjacency(self):
        omega_0 = 2*math.pi*1e3
        for i in range(self.N):
            if i == 0:
                self.matrix[i, 1] = np.exp(1j * math.pi*np.sin(self.teta[0]))
                self.matrix[i, -1] = np.exp(1j * math.pi*np.sin(self.teta[0]) * (self.N - 1))
            elif i == self.N - 1:
                self.matrix[i, -2] = np.exp(-1j * math.pi*np.sin(self.teta[0]))
                self.matrix[i, 0] = np.exp(-1j * math.pi*np.sin(self.teta[0]) * (self.N - 1))
            else:
                self.matrix[i, i - 1] = np.exp(-1j * math.pi*np.sin(self.teta[0]))
                self.matrix[i, i + 1] = np.exp(1j * math.pi*np.sin(self.teta[0]))
        A_time = 0.5 * self.matrix

        for i in range(self.K):
            if i == 0:
                self.matrix2[i, 1] = np.exp(-1j * omega_0)
                self.matrix2[i, -1] = np.exp(-1j * omega_0 * (self.K - 1))
            elif i == self.K - 1:
                self.matrix2[i, -2] = np.exp(1j * omega_0)
                self.matrix2[i, 0] = np.exp(1j * omega_0 * (self.K - 1))
            else:
                self.matrix2[i, i - 1] = np.exp(1j * omega_0)
                self.matrix2[i, i + 1] = np.exp(-1j * omega_0)
        A_space = 0.5 * self.matrix2

        return np.kron(A_space, A_time)

if __name__ == "__main__":
    print("Not main file")