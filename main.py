from classes import prameters_class, Matrix_class
from functions import observ, GFT, quantize, G_DOA,angles_generate,G_DOA2
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__": #TODO RMSE (monte)
    N_a = [10]
    N_q = [0]
    D = 2
    teta_range = [-60, 60]
    # SNR = 0
    SNR_space = np.linspace(-5, 5, 5)
    snap = 100
    # snap_space = np.linspace(100, 1000, 10, dtype=int)
    monte = 10
    delta = 5 #Minimal gap between two determenistic angles
    Res = 1
    # delta_space = np.linspace(0.8, 6, 20)
    relevant_space = SNR_space  # TODO
    Error1 = np.zeros((len(relevant_space), len(N_q)))
    for i in range(len(relevant_space)):
        for j in range(len(N_q)):
            my_parameters = prameters_class(M=N_a[j] + N_q[j],N_q=N_q[j], SNR=SNR_space[i], K=snap, D=D,teta_range=teta_range, monte=monte,
                                            delta=delta,Res=Res)
            Error1[i, j] = G_DOA(my_parameters)
    np.save("Error1.npy", Error1)
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


