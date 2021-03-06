Mesh Fracturing
===============

The following is a basic mesh fracturing application built in openGL. The program leverages the use of a geometry shader in order to create individual fragments of an object. The geometry shader calculates the normals of each individual fragment and translates the fragment in the direction of its normal.

Building The Project
--------------------
The project can be built using the supplied CMakeLists file. The src folder contains all relevant C++ code. The resources folder, which the program takes as a single argument, contains all obj files, as well as shaders that the program utilizes.

Technologies Used
-----------------
* Texture Mapping
* Pitch & Yaw
* Collision Detection
* Geometry Shader

Features
--------
* Randomly generated world of 100 objects
* Free camera movement
* Mouse click detection using axis aligned bounding boxes
* User can control time of the world

Controls
--------
* WASD and mouse to move around the world
* Click on a mesh to fracture it
* Right and Left arrow keys control "time"
* Up and Down arrow keys control the speed of "time"
* Spacebar stops "time"
