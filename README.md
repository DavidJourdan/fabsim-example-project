# Fabsim example project

Uses the optim library for static solves and Polyscope for the viewer

## Build instructions

``` 
git clone --recursive https://github.com/DavidJourdan/fabsim-example-project 
cd fabsim-example-project 
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Then you can try ```cd build; ./rod_example``` for an example using the discrete elastic rods model or ```cd build; ./membrane_example``` for a membrane simulation using the StVK FEM model (other material models are available as well)
