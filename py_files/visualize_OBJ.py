import igl
import meshplot as mp
import numpy as np


#mp.website() # plot in 3d viewer using Html 
mp.Viewer.to_html() 

# loading and plotting mesh
v, f = igl.read_triangle_mesh('./data/processed_data/0002_letti.obj')
k = igl.gaussian_curvature(v, f)
plot_mesh = mp.plot(v, f, return_plot = True)
plot_mesh1 = mp.plot(v,f, k, return_plot=True)
print('shape of vertices using igl library = ',v.shape)
print('shape of triangles using igl library = ', f.shape)

# d = mp.subplot(v, f, c=v[:, 1], s=[2, 2, 0])
# mp.subplot(v, f, c=n, s=[2, 2, 1], data=d)
# mp.subplot(v, f, c=np.random.rand(*f.shape), s=[2, 2, 2], data=d)
# mp.subplot(v, f, c=fs, s=[2, 2, 3], data=d)