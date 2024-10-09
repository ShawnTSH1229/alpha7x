# Alpha7XRender

Alpha 7X is a tiny offline renderer for graphics learning. It references PBRT and uses embree to accelerate the ray-scene intersection. Alpha7x support traditional path tracing and stochastic progressive photon mapping.

# Getting Started

Visual Studio 2019 or 2022 is recommended, XEngine is only tested on Windows.

1.Cloning the repository with `git clone https://github.com/ShawnTSH1229/alpha7x.git`.

2.Configuring the build

```shell
# Create a build directory
mkdir build
cd build

# x86-64 using a Visual Studio solution
cmake -G "Visual Studio 17 2022" ../
```
# Example

stochastic progressive photon mapping ([my related post link](https://shawntsh1229.github.io/2024/10/06/Stochastic-Progressive-Photon-Mapping-In-Alpha7XRenderer/)):

<p align="left">
    <img src="/resource/result.png" width="50%" height="50%">
</p>