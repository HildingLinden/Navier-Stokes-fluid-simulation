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

## C++ Multithreaded AVX2
This version is implemented in c++. It uses c++11 threads for multithreading and avx2 intrinsics for vectorization.  The program starts a pool of threads at the start of the program. Each thread then waits on a condition_variable until work is available. The work is a std::function lambda that calls advectLoop, projectHeightMapLoop, projectMassConservLoop or linearSolveLoop with the indices of the loop that the thread should run. When the thread is done it notifies the main thread and goes back to waiting for work. When the main thread is notified and all threads have completed it will continue the program. Storing the result in the array does not need to be synchronized since each thread writes to it's own region of memory. However, the performance can suffer a bit if two threads write to regions within the same chache line (within a couple of elements from each other) at the same time. This is very unlikely to happen and it is only a couple of elements per thread that is affected so the hit to the total performance is negligable.
6 threads seem to perform best on my 6700k quad-core processor when the resolution of the fluidgrid is atleast 200x200.

## CUDA
This version will be implemented in CUDA to be able to run on NVIDIA GPUs.

## Vulkan
This version will be implemented as a Vulkan compute shader and will also be using Vulkan to display the results.