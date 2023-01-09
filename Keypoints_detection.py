import igl
import meshplot as mp
import numpy as np

# loading and plotting mesh
vert, tri = igl.read_triangle_mesh('./data/processed_data/0002_letti.obj')

# to compute FPFH feature descriptor to detect keypoints
def compute_fpfh(mesh, support_radius):
  fpfh = []
  for i, (point, normal) in enumerate(mesh):
    histogram = np.zeros(NUM_BINS) # initialize histogram to zero
    # find neighboring points within support radius
    neighbors = []
    for j, (neighbor, _) in enumerate(mesh):
      if np.linalg.norm(point - neighbor) < support_radius:
        neighbors.append(j)
    # calculate histogram
    for j in neighbors:
      neighbor, neighbor_normal = mesh[j]
      # calculate angle between normals
      angle = np.arccos(np.dot(normal, neighbor_normal))
      # increment appropriate histogram bin
      bin_idx = int(angle / BIN_SIZE)
      histogram[bin_idx] += 1
    # normalize histogram
    histogram /= len(neighbors)
    fpfh.append(histogram)
  return fpfh
