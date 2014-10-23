HE-tumor-object-segmentation
============================

Segmentation of tumor objects in H&amp;E images using graph-based clustering based on opponent color space and internal image statistics
Author: Luong Nguyen
Email: lun5@pitt.edu
Started 10/23/2014

Main components:
1. Inputs: tiling of WSI, parallelization, rotation matrix
2. Annotations: collect annotations from pathologists
3. Affinity calculation: Features and affinity function
4. Clustering: NMinCut, Embedded Spectral Clustering 

Reused codes and packages from 
Citation of KDE
Citation of Multiscale Spectral Clustering 
Citation of Crisp Boundary 

