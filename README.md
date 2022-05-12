# SurfaceFlowSimulatorCuda
Project for accurately and rapidly simulating the surface flow discharge

# Run-time environment
Visual Studio and CUDA8.0 are installed on a computer with the NVIDIA graphics card which the compute capability is 2.1 or above

#detail procedure of the parallel project

1)Downloading all of the files except the ReadMe.txt from the website (https://github.com/parallelProject/SurfaceFlowSimulatorCuda)
  (Tip:the gdal.lib file in the "lib" folder of the "include" folder is larger than 50M, it must be downloaded in seperate and replace the gdal.lib file in the *.ZIP file.)
 
2)Setting the file directory of all the variables in the SurfaceFlowSimulator.cpp and the pre-defined parameters
  
3)Modifing the variables "DOWNLOADFILE_PATH"  to the your path

4)Running  and finally obtaining the 365 images with the daily surface flow discharge and a txt file recording the computing time
