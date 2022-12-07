import igl
import meshplot as mp
import numpy as np



# loading and plotting mesh
vert, tri = igl.read_triangle_mesh('./data/processed_data/0002_letti.obj')
#plot_mesh = mp.plot(vert, tri)
print('shape of vertices = ',vert.shape)
print('shape of triangles', tri.shape)