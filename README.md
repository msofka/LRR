
### Introduction

The code implements:

Sofka, M., V. Stewart, C., 2010. [Location Registration and Recognition (LRR) for Serial Analysis of Nodules in Lung CT Scans](http://www.sofka.com/pdfs/sofka-mia10.pdf). Medical Image Analysis 14, 407–428.

<p align="center"><img align="center" src="https://github.com/msofka/LRR/blob/master/LRR_aligned.png" width="500"/></p>
Figure: A panel of nine images for each result presented to an independent observer for the alignment evaluation. The rows show axial, coronal, and sagittal slices. The columns show mapped moving slices, fixed slices, and checkerboards alternating mapped moving and fixed slices. The features are superimposed onto the images.


### Build

CMake build system is recommended for creating make files and project files. CMake files are included with the source code distribution.

### Run

The algorithm components are run separately without including them in one executable. This makes the components easier to debug, evaluate, and integrate in other projects. The disadvantage is that data is loaded and saved multiple times. What follows are the steps to run the LRR algorithm. The steps consists of calling each executable from simple bash scripts (use directly in Linux or download Cygwin to use on Windows). In the scripts, you will need to modify path to your binaries. You can also run the executables directly (see below). 

# Step 1: Feature Extraction

Download example volume [here](http://www.cs.rpi.edu/~sofka/cgi-bin/downloadvw.cgi?volume_pair).

Extract features, keypoints, and descriptors:

```
 ./extract_features.bash volumes_with_nodules.txt
```

The script also calls a script that computes a Voronoi map for each feature set:

```
./compute_distance_map.bash volumes_with_nodules.txt
```

In addition, the script also computes watershed oversegmentations by running:

```
 ./watershed_segmentation.bash volumes_with_nodules.txt
```

Features, keypoints, and descriptors are stored in a VTK format.

Alternatively, you can run the command line executables directly for each volume:

```
 cd 0009/4893
 ExtractFeatures 4978whole 4978whole_00.vtk
 ExtractKeypoints 4978whole_00.vtk 4978wholekeypoints.vtk
 ComputeDescriptors 4978wholekeypoints.vtk 4978whole_00.vtk 4978wholedesc.vtk
 compute_distance_map 4978whole
 WatershedSegmentation1 4978whole 4978wholewatershed 2.0 10 0.001 0.10
```

Repeat the steps for volume in 0009/4717. 

# Step 2: Initialization
Generate initializations by keypoint indexing:

```
 ./invariant_indexing.bash pairs_with_nodules.txt
```

This creates a directory `nodules` which contains: 1) moving and fixed descriptors for each match, 2) images surrounding fixed and moving keypoints, 3) the query location and its neighborhood within fixed image. 

Running the executable directly:

```
 cd 0009
 indexing_one_descriptor.exe 4717/4722wholedesc.vtk 4893/4978wholedesc.vtk \
  -mov 4717/4722whole -fix 4893/4978whole \
  -trans 4893_4978whole-4717_4722wholediff.vtk -locs  4893/4978nodules.txt
```

(Note: Images are specified for debugging purposes. They are not used in the algorithm.) 

# Step 3: Estimation
Use the initializations generated in the previous step to refine the estimate and output the verified transform:

```
 ./location_registration.bash pairs_with_nodules.txt
```

This creates a directory ``nodulesreg'' that contains: 1) a text file with the final estimated transform and various statistics, and 2) a transform and its inverse in VTK format. If LAST_ITER is defined as 1 in the code, the executable will also produce a panel of 9 images containing fixed, moving, and checkerboard images of axial, sagittal, and coronal slices. 
Running the executable directly:

```
 cd 0009     
 location_registration.exe 4893/4978whole 4717/4722whole \
 -trans 4893_4978whole-4717_4722wholediff.vtk \
 -matches 4893/4978nodules.txt -segment watershed
 ```
(Note: Images are specified for debugging purposes. They are not used in the algorithm.)
