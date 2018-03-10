# CS130-pipeline

In this assignment, I will be implementing a simplified 3D rendering pipeline (with flat shading). This will consist of several parts, introduced through a series of tests:

  1. Vertex and viewing transformations
  2. Rasterization and interpolation
  3. Clipping
  4. Using a z-buffer for hidden surfaces
  
## How to run it

Download, unzip, and do make file first (enter make in terminal)

Then if you're running this over SSH, without the X server, you can use:

  `./minigl-nogl -i tests/xx.txt -s tests/xx.png -p`

Otherwise, you can try this:

  `./minigl -g -i tests/xx.txt`

If you just wanna see the grading, you can use:

  `chmod -x grading-helper.pl` and then `bash grading-script.sh tests`
