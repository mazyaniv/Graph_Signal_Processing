import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
import scipy.signal as ss
from classes import prameters_class, Matrix_class
from functions import observ, GFT, quantize

if __name__ == "__main__":
    theta = [-15,32]
    my_parameters = prameters_class(20,10,200,theta)
    steering_original = Matrix_class(my_parameters).steering()
    obs_a = observ(my_parameters.SNR, my_parameters.K, steering_original)
    x_vec = quantize(obs_a, my_parameters.M)

    #Method 1:
    # x_vec_s = x_vec.T.flatten()
    # A_s_original = Matrix_class(my_parameters).Adjacency()
    # x_tag_s_original = GFT(A_s_original,x_vec_s) #[np.abs(eigenvectors[:,i].T.conjugate()@np.tile(np.exp(-1j * np.pi * np.arange(4) * np.sin(np.radians(theta))),3)) for i in range(12)]
    #
    # my_parameters2 = prameters_class(my_parameters.M, my_parameters.SNR, my_parameters.K, my_parameters.D, [1])
    # A_s = Matrix_class(my_parameters2).Adjacency()
    # x_tag_s = GFT(A_s,x_vec_s)
    #
    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s_original) / np.max(np.abs(x_tag_s_original)), use_line_collection=True) #look only at \hat{x} abs
    # plt.show()
    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s)/np.max(np.abs(x_tag_s)), use_line_collection=True)
    # plt.show()

    # print(1/LA.norm(np.abs(np.delete(x_tag_s_original, np.argmax(np.abs(x_tag_s_original))))/np.max(np.abs(x_tag_s_original))))
    # print(1 / LA.norm(np.abs(np.delete(x_tag_s, np.argmax(np.abs(x_tag_s)))) / np.max(np.abs(x_tag_s))))

    # Method 2:
    # x_vec_s = np.average(x_vec,1)
    # A_s_original = steering_original @ steering_original.T.conjugate() - (1/(my_parameters.M - 1))*np.eye(my_parameters.M)
    # x_tag_s_original = GFT(A_s_original, x_vec_s)
    #
    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s_original) / np.max(np.abs(x_tag_s_original)), use_line_collection=True) #look only at \hat{x} abs
    # plt.show()
    # my_parameters2 = prameters_class(my_parameters.M, my_parameters.SNR, my_parameters.K,[70])
    # steering = Matrix_class(my_parameters2).steering()
    # A_s = steering @ steering.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)
    # x_tag_s = GFT(A_s, x_vec_s)
    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s)/np.max(np.abs(x_tag_s)), use_line_collection=True)
    # plt.show()

    def G_DOA(pram,teta_range):
        x_vec_s = np.average(x_vec,1)#x_vec.T.flatten()
        theta_range = np.arange(teta_range[0], teta_range[1], pram.Res)
        spectrum = np.zeros(len(theta_range))
        for idx, theta in enumerate(theta_range):
            my_parameters2 = prameters_class(pram.M, pram.SNR, pram.K, [theta])
            steering2 = Matrix_class(my_parameters2).steering()
            A_s = steering2 @ steering2.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)#Matrix_class(my_parameters2).Adjacency()
            x_tag_s = GFT(A_s, x_vec_s)
            spectrum[idx]= 1/LA.norm(np.abs(np.delete(x_tag_s, np.argmax(np.abs(x_tag_s))))/np.max(np.abs(x_tag_s)))
        plt.plot(theta_range, spectrum,marker="*")
        plt.show()

        peaks, _ = ss.find_peaks(spectrum)
        peaks = list(peaks)
        peaks.sort(key=lambda x: spectrum[x])
        pred = np.array(peaks[-pram.D:])
        pred = np.sort(pred)[::-1]
        return pred*pram.Res+teta_range[0]

teta_range = [-60, 60]
print(G_DOA(my_parameters,teta_range))






