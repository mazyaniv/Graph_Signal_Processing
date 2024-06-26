import numpy as np
import math
class prameters_class():
    def __init__(self, M,N_q,SNR,K,D,teta_range,monte=100,delta=5,Res=1):
        self.M = M #M
        self.N_q = N_q
        self.SNR = SNR
        self.K = K
        self.D = D
        self.teta_range = teta_range
        self.monte = monte
        self.delta = delta
        self.Res = Res
class Matrix_class():
    def __init__(self, pram):
        self.M = pram.M #sensors
        self.K = pram.K #snapshots
        self.D = pram.D #sources
        self.Stee = np.zeros((self.M, self.D), dtype=complex) #steering matrix
        self.matrix = np.zeros((self.M*self.K, self.M*self.K), dtype=complex)
    def steering(self,theta):
        A_mask = np.zeros((self.M, self.D), dtype=complex)
        for j in range(self.D):
            A_mask[:, j] = np.exp(-1j * np.pi * np.arange(self.M) * np.sin(np.radians(theta[j])))
        self.Stee = A_mask
        return self.Stee
    # def Adjacency(self):
    #     matrix1 = np.zeros((self.M, self.M), dtype=complex)
    #     matrix2 = np.zeros((self.K, self.K), dtype=complex)
    #     omega_0 = 2*math.pi*1e3
    #     for i in range(self.M):
    #         if i == 0:
    #             matrix1[i, 1] = np.exp(1j * math.pi*np.sin(np.radians(self.teta[0])))
    #             matrix1[i, -1] = np.exp(1j * math.pi*np.sin(np.radians(self.teta[0])) * (self.M - 1))
    #         elif i == self.M - 1:
    #             matrix1[i, -2] = np.exp(-1j * math.pi*np.sin(np.radians(self.teta[0])))
    #             matrix1[i, 0] = np.exp(-1j * math.pi*np.sin(np.radians(self.teta[0])) * (self.M - 1))
    #         else:
    #             matrix1[i, i - 1] = np.exp(-1j * math.pi*np.sin(np.radians(self.teta[0])))
    #             matrix1[i, i + 1] = np.exp(1j * math.pi*np.sin(np.radians(self.teta[0])))
    #     A_time = 0.5 * matrix1
    #
    #     for i in range(self.K):
    #         if i == 0:
    #             matrix2[i, 1] = np.exp(-1j * omega_0)
    #             matrix2[i, -1] = np.exp(-1j * omega_0 * (self.K - 1))
    #         elif i == self.K - 1:
    #             matrix2[i, -2] = np.exp(1j * omega_0)
    #             matrix2[i, 0] = np.exp(1j * omega_0 * (self.K - 1))
    #         else:
    #             matrix2[i, i - 1] = np.exp(1j * omega_0)
    #             matrix2[i, i + 1] = np.exp(-1j * omega_0)
    #     A_space = 0.5 * matrix2
    #     self.matrix = np.kron(A_time,A_space)
    #     return self.matrix

if __name__ == "__main__":
    print("Not main file")