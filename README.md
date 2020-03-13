# Navier-Stokes fluid simulation
In this repository I will implement a fluid simulator + visualisation based on the paper [Real-Time Fluid Dynamics for Games](https://pdfs.semanticscholar.org/847f/819a4ea14bd789aca8bc88e85e906cfc657c.pdf) by Jos Stam. There will be multiple implementations using different programming languages and APIs.  
This project is mostly for learning the basics in CFD and optimization of computationally heavy code.

More resources:
https://mikeash.com/pyblog/fluid-simulation-for-dummies.html  
http://developer.download.nvidia.com/books/HTML/gpugems/gpugems_ch38.html  
https://developers.google.com/web/updates/2019/08/get-started-with-gpu-compute-on-the-web

## Javascript only
This version is implemented in Javascript only and uses the Canvas 2D context for rendering the result.  
[Demo here](https://hildinglinden.github.io/Navier-Stokes-fluid-simulation/javascript/)

## WebGPU
This version will also be implemented in Javascript but will also use the WebGPU API which will allow the code to run on the GPU.  
[WebGPU implementation status](https://github.com/gpuweb/gpuweb/wiki/Implementation-Status)

## C/C++
This version will be implemented in standard C or C++.

## Vectorized and multithreaded C
This version will also be implemented in C but with AVX2 vectorization and multithreading with OpenMP.

## CUDA
This version will be implemented in CUDA to be able to run on NVIDIA GPUs.