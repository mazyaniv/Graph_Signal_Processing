from classes import prameters_class, Matrix_class
from functions import observ, GFT, quantize, G_DOA,angles_generate
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__": #TODO RMSE (monte)
    my_parameters = prameters_class(M=20,SNR=-5,K=100,teta_range=[-60, 60],D=2,monte=10,Res=1)
    # angles = angles_generate(my_parameters)
    # print(angles)
    # S = Matrix_class(my_parameters).steering(angles)
    print(G_DOA(my_parameters,0))

