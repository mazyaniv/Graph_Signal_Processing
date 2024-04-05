from classes import prameters_class, Matrix_class
from functions import observ, GFT, quantize, G_DOA,angles_generate,G_DOA2
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
    # my_parameters = prameters_class(M=20, N_q=0, SNR=0, K=1000, D=2, teta_range=[-60, 60],monte=1,delta=5,Res=0.5)
    # print(G_DOA(my_parameters))
    # theta = [20,30]
    # steering_original = Matrix_class(my_parameters).steering(theta)
    # obs_a = observ(my_parameters.SNR, my_parameters.K, steering_original)
    # x_vec = quantize(obs_a, 0)
    # x_vec_s = np.average(x_vec,1)
    #
    # steering = np.exp(-1j * np.pi * np.arange(my_parameters.M) * np.sin(np.radians(theta[0]))).reshape(my_parameters.M,1)
    # A_s_original = steering_original @ steering_original.T.conjugate() - (1/(my_parameters.M - 1))*np.eye(my_parameters.M)
    # x_tag_s_original = GFT(A_s_original, x_vec_s)
    #
    # theta_fake = [-10]
    # steering2 = np.exp(-1j * np.pi * np.arange(my_parameters.M) * np.sin(np.radians(theta_fake))).reshape(my_parameters.M,1)
    # A_s = steering2 @ steering2.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)
    # x_tag_s = GFT(A_s, x_vec_s)
    #
    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s_original),linefmt='b-',label="correct")
    # plt.stem(np.abs(x_tag_s),linefmt='g-',markerfmt='gx',label="Incorrect")
    # # plt.title("Incorrect")
    # plt.legend()
    # plt.show()
###############################################################
    N_a = [20, 2, 0]
    N_q = [0, 18, 20]
    D = 2
    teta_range = [-60, 60]
    # SNR = 0
    SNR_space = np.linspace(-5, 5, 5)
    snap = 1000
    # snap_space = np.linspace(100, 1000, 10, dtype=int)
    monte = 100
    delta = 5 #Minimal gap between two determenistic angles
    Res = 0.5
    # delta_space = np.linspace(0.8, 6, 20)
    relevant_space = SNR_space  # TODO
    Error1 = np.zeros((len(relevant_space), len(N_q)))
    for i in range(len(relevant_space)):
        for j in range(len(N_q)):
            my_parameters = prameters_class(M=N_a[j] + N_q[j],N_q=N_q[j], SNR=SNR_space[i], K=snap, D=D,teta_range=teta_range, monte=monte,
                                            delta=delta,Res=Res)
            Error1[i, j] = G_DOA(my_parameters)
    # print(Error1)
    np.save(f"RMSE for delta={my_parameters.delta}, Monte={my_parameters.monte}, Res={my_parameters.Res},Snap={my_parameters.K}.npy", Error1)
    fig = plt.figure(figsize=(12, 8))
    colors = ['red', 'b', 'black']
    for i in range(len(N_q)):  # TODO
        plt.plot(relevant_space, Error1[:, i], linestyle='solid', marker=".", color=colors[i],
                 label=f'N_a={N_a[i]},N_q={N_q[i]}')
    plt.grid()
    plt.title(f"RMSE for $\Delta$={my_parameters.delta}, Monte={my_parameters.monte}, Res={my_parameters.Res}, "
              f"Snap={my_parameters.K}")  # TODO
    plt.ylabel("RMSE")
    plt.xlabel(r"$SNR_{dB}$")
    plt.legend(loc='upper right', fontsize='small')
    plt.show()


