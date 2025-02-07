# InnoSaddle

### Overview

InnoSaddle is a project focused on utilizing depth data from Apple iOS true depth sensors to detect key points in depth maps and perform measurements at specific positions. This repository contains the code and documentation for my master thesis, which explores the application of 3D machine learning techniques in processing RGBD data.

### Project Objectives

- **Depth Data Acquisition**: Collect RGBD data using Apple iOS true depth sensors.
- **Key Point Detection**: Implement algorithms for detecting key points in depth maps.
- **Measurement Analysis**: Develop methods for measuring distances or other metrics at required positions based on detected key points.

### Libraries and Tools Used

- **Open3D**: For 3D data processing and visualization.
- **PyTorch3D**: For 3D deep learning tasks, including mesh manipulation and rendering.
- **Panda3D**: For 3D visualization and simulation (optional).

### Steps to Replicate

1. **Data Collection**:
   - Use Apple iOS devices with true depth sensors to capture RGBD data.
   - Store the data in a suitable format for processing.

2. **Key Point Detection**:
   - Implement key point detection algorithms on the depth maps.
   - Use libraries like Open3D for efficient processing.

3. **Measurement Analysis**:
   - Develop scripts to measure distances or other metrics based on detected key points.
   - Utilize PyTorch3D for advanced 3D analysis if necessary.

4. **Visualization (Optional)**:
   - Use Mayavi or pi3d for interactive 3D visualization of results.

5. **Conversion of 3D Models**:
   - Convert USDZ files to OBJ or PLY formats.
   - Use Open3D to read and manipulate these files.

### Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/surya-1729/InnoSaddle.git
   ```

2. Install required libraries:
   ```bash
   pip install open3d torch torchvision
   ```

3. Run key point detection and measurement scripts:
   ```bash
   python detect_keypoints.py
   python measure_positions.py
   ```
