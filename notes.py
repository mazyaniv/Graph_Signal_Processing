import numpy as np



# obs_a = observ(my_parameters.SNR, my_parameters.K, steering_original)
# x_vec = quantize(obs_a, 0)
# Method 1:
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
#
# print(1/LA.norm(np.abs(np.delete(x_tag_s_original, np.argmax(np.abs(x_tag_s_original))))/np.max(np.abs(x_tag_s_original))))
# print(1 / LA.norm(np.abs(np.delete(x_tag_s, np.argmax(np.abs(x_tag_s)))) / np.max(np.abs(x_tag_s))))

# Method 2:
# x_vec_s = np.average(x_vec,1)
# A_s_original = steering_original @ steering_original.T.conjugate() - (1/(my_parameters.M - 1))*np.eye(my_parameters.M)
# x_tag_s_original = GFT(A_s_original, x_vec_s)
#
# my_parameters2 = prameters_class(my_parameters.M, my_parameters.SNR, my_parameters.K,[70])
# steering = Matrix_class(my_parameters2).steering()
# A_s = steering @ steering.T.conjugate() - (1 / (my_parameters.M - 1)) * np.eye(my_parameters.M)
# x_tag_s = GFT(A_s, x_vec_s)
#
# # plt.figure(figsize=(10, 6))
# # plt.stem(np.abs(x_tag_s_original) / np.max(np.abs(x_tag_s_original)), use_line_collection=True)  # look only at \hat{x} abs
# # plt.title("Original")
# # plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.stem(np.abs(x_tag_s_original),linefmt='b-',label="correct")
# plt.stem(np.abs(x_tag_s),linefmt='g-',markerfmt='gx',label="Incorrect")
# # plt.title("Incorrect")
# plt.legend()
# plt.show()