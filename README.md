# Fabsim example project

Uses the optim library for static solves and Polyscope for the viewer

## Build instructions

``` 
git pull --recursive https://github.com/DavidJourdan/fabsim-example-project 
cd fabsim-example-project 
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Then you can try ```build/rod_example``` or ```build/membrane_example```
