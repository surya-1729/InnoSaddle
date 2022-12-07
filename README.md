# InnoSaddle

Depth data from true depth sensor:

We get the RGBD data from Apple IOS true depth sensor. Now, have to look for key point detection in the depth maps and make the measurements for our required positions.

Libraries for 3D Machine Learning:

1. Panda3D

    An open-source library specifically used for 3D games, simulations, and visualizations. Written in C++ along with Python bindings, this library is accountable for gloss mapping, normal mapping, cartoon shading and inking, and HDR, etc.

    Best features:

    A cross-platform engine that eases deployment on most of the supported platforms.
    It also combines command-line tools that help in optimizing and processing source assets.
    In combination with C++ and Python, you get a faster rate of development without needing to compromising on the performance.

2. PyTorch3D

    Another open-source library used for 3D deep learning. It is highly optimized and highly modular and is specially designed to make 3D deep learning much easier using the PyTorch library. This library is used by Facebook AI research to boost research projects like Mesh R-CNN.

    Best features:

    Provides data structure to help store and manipulate triangle meshes.
    PyTorch3D offers a differentiable mesh renderer.
    Supported by a whole bunch of cloud platforms that provide easy scaling and frictionless development.


3. Open3D

    An open-source library ideally used for 3D data processing which also supports the rapid development of software that deals with 3D data. The front-end of Open3D is exposed to a set of selected data structures and algorithms in both the languages – C++ and Python. However, the back-end is highly optimized and the set-up is used for parallelization.

    Best features:

    3D visualization
    3D data processing algorithms
    3D data structures
    Physically-based rendering (PBR)
    Scene reconstruction

4. Mayavi

    A cross-platform tool used for 3D scientific data visualization. It is open-source and is completely written in Python language. Mayavi provides a clean scripting interface to have an easy and interactive visualization of 3D data. Besides a clean scripting interface, this library also provides an object-oriented programming interface, one-liners, and many more other features.

    Best features:

    Offers easy extendibility through modules, custom sources, and data filters.
    Helps in reading multiple file-formats including PLOT3D and VTK (legacy and XML).
    Saves rendered visualizations in multiple image formats.

5. pi3d

    It is a Python module that aims at simplifying writing 3D in Python language while also allowing them access to raspberry Pi GPU. The pi3d easily enables both 2D and 3D rendering while providing access to host multiple commands to load animated models and textured models, and to create shaders and fractal landscapes, etc.

    Best features:

    Apart from Raspberry Pi GPU, pi3d also runs well on Linux using X server or Windows using Pygame and on Android using python for android.
    pi3d offers compatibility features in both Python 2 and Python 3.

STEPS further:

1. Convert :USDZ to .OBJ or .PLY
2. Use this .OBJ or .PLY file to read points information using open3d
3. 