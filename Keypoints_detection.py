import open3d as o3d

source1 = o3d.io.read_triangle_mesh('./data/processed_data/0002_letti.obj')
print(source1)

source2 = o3d.io.read_point_cloud('./data/processed_data/0002_letti.obj')
print(source2)