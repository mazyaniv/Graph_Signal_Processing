from classes import prameters_class, Matrix_class
from functions import observ, GFT, quantize, G_DOA,angles_generate,G_DOA2, general
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
    # my_parameters = prameters_class(M=20, N_q=0, SNR=20, K=500, D=2, teta_range=[-100, 100],monte=1,delta=20,Res=1)
    # print(G_DOA(my_parameters))
    # theta = [-30,45]
    # steering_original = Matrix_class(my_parameters).steering(theta)
    # obs_a = observ(my_parameters.SNR, my_parameters.K, steering_original)
    # x_vec = quantize(obs_a, 0)
    # x_vec_s = np.average(x_vec,1)
    #
    # steering1 = np.exp(-1j * np.pi * np.arange(my_parameters.M) * np.sin(np.radians(theta[0]))).reshape(my_parameters.M,1)
    # A_s_original1 = steering1 @ steering1.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)
    # x_tag_s_original = GFT(A_s_original1, x_vec_s)
    #
    # steering2 = np.exp(-1j * np.pi * np.arange(my_parameters.M) * np.sin(np.radians(theta[1]))).reshape(my_parameters.M,1)
    # A_s_original2 = steering2 @ steering2.T.conjugate() - (1/(my_parameters.M - 1))*np.eye(my_parameters.M)
    # x_tag_s_original2 = GFT(A_s_original2, x_vec_s)
    #
    # theta_fake = [10]
    # steering2 = np.exp(-1j * np.pi * np.arange(my_parameters.M) * np.sin(np.radians(theta_fake))).reshape(my_parameters.M,1)
    # A_s = steering2 @ steering2.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)
    # x_tag_s = GFT(A_s, x_vec_s)

    # plt.figure(figsize=(10, 6))
    # plt.stem(np.abs(x_tag_s_original),linefmt='b-',label=r"$\theta={t}^\degree$".format(t=theta[0]))
    # plt.stem(np.abs(x_tag_s_original2),linefmt='r--',markerfmt='rx',label=r"$\theta={t}^\degree$".format(t=theta[1]))
    # plt.stem(np.abs(x_tag_s),linefmt='g-',markerfmt='g.',label=r"$\theta={t}^\degree$".format(t=theta_fake[0]))
    # plt.xlabel("GFT eigenvalue index")
    # plt.ylabel("GFT cofficient magnitude")
    # plt.title(r"DOA for $\theta={theta}^\degree$".format(theta=theta))
    # plt.legend()
    # plt.show()
###############################################################
    N_a = [20]
    N_q = [0]
    D = 2
    teta_range = [-60, 60]
    # SNR = 0
    SNR_space = np.linspace(-15, -5, 8)
    snap = 153
    # snap_space = np.linspace(100, 1000, 10, dtype=int)
    monte = 500
    delta = 20 #Minimal gap between two determenijjstic angles
    Res = 1
    # delta_space = np.linspace(0.8, 6, 20)
    relevant_space = SNR_space  # TODO

    Error1 = np.zeros((len(relevant_space), len(N_q)))
    Error2 = np.zeros((len(relevant_space), len(N_q)))
    for i in range(len(relevant_space)):
        for j in range(len(N_q)):
            my_parameters = prameters_class(M=N_a[j] + N_q[j],N_q=N_q[j], SNR=SNR_space[i], K=snap, D=D,teta_range=teta_range, monte=monte,
                                            delta=delta,Res=Res)
            Error1[i, j] = G_DOA(my_parameters)
            Error2[i, j] = general(my_parameters)
    print("RMSE=",Error1)
    # np.save(f"RMSE for Monte={my_parameters.monte},Snap={my_parameters.K}.npy", Error1)
    # Error1 = np.load(f"RMSE for Monte={monte},Snap={snap}.npy")
    fig = plt.figure(figsize=(12, 8))
    colors = ['black', 'red']
    for i in range(len(N_q)):  # TODO
        plt.plot(relevant_space[2:-2], Error1[2:-2, i], linestyle='solid', marker=".", color=colors[i])#,
                 # label=f'N_a={N_a[i]},N_q={N_q[i]}')
        # plt.plot(relevant_space, Error2[:, i], linestyle='dashed', marker="x", color=colors[i],
        #          label=f'N_a={N_a[i]},N_q={N_q[i]},MUSIC')
    ax = plt.gca()
    ax.set_xticks(np.arange(SNR_space[2], SNR_space[-2], 0.25), minor=True)
    ax.grid(which='major', alpha=0.5)
    ax.grid(which='minor', linestyle="--", alpha=0.25)
    plt.title(f"RMSE for {monte} simulations and {snap} snapshots")  # TODO
    plt.ylabel("RMSE")
    plt.xlabel(r"$SNR_{dB}$")
    plt.xlim(SNR_space[2], SNR_space[-3])
    # plt.legend(loc='upper right', fontsize='small')
    plt.show()

