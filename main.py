import numpy as np
from matplotlib import pyplot as plt
from classes import prameters_class, Matrix_class
from functions import observ, GFT

if __name__ == "__main__":

    D = 1 #TODO
    N = 3 #sensors
    teta = [20]
    SNR = 0
    k = 4 #samples

    monte = 100
    delta = 5 #Minimal gap between two determenistic angles
    Res = 1

    my_parameters = prameters_class(N,SNR,k,D,teta,monte,delta,Res)
    matrix = Matrix_class(my_parameters.N, my_parameters.k, teta)
    A = matrix.Adjacency()
    steering = matrix.steering()
    x_vec = observ(my_parameters.SNR,my_parameters.k,steering)
    x_tag = GFT(A,x_vec)

    plt.figure(figsize=(10, 6))
    plt.stem(np.abs(x_tag), use_line_collection=True)
    plt.show()



