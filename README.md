# Planar triangulation for N-view

Manuscript in review

---
This repository contains the code 
for the N-view planar triangulation 
solver explained in [this paper](NONE) [1]. 


**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GPLv3](https://github.com/mergarsal/FastNViewTriangulation/blob/main/LICENSE)


If you use this code for your research, please cite:

```
NONE
```

## Dependencies

- Eigen 



## Build

```
git clone https://github.com/mergarsal/NPlanarCertTriangulation.git
cd NPlanarCertTriangulation

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. Run: 

```
        ./bin/example_base
```
 


## Install 
In `build` folder: 

```
        sudo make install
```

We also provide the uninstall script: 

```
        sudo make uninstall
```





# Use in external project 
1. Install our library with 

```
sudo make install 
```

2. In your project, find the library. 


```
find_package(NPlanarTrian REQUIRED)
```

3. Include the library as target, e.g., 

```
add_executable(example_base ${CMAKE_CURRENT_SOURCE_DIR}/example_base.cpp)

target_link_libraries(example_base 
                              NPlanarTrian      
                     )
```              

