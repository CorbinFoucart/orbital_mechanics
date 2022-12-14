# orbital_mechanics

## time integration of the n-body problem.

![](https://github.com/CorbinFoucart/orbital_mechanics/blob/main/three_body.gif)

- Simulation code is written in C++, see `celestial_mechanics.h`,
`celestial_mechanics.cc`, and the driver code `driver.cc`. Outputs trajectories
to CSV file.
- Simulation code is templated by space dimension to simultaneously implement
  both 2D and 3D dynamics in a dimension-independent way.
- Python code to plot and animate the CSV trajectory output can be found in
  `animations.ipynb`.

[Some animations.](https://www.youtube.com/playlist?list=PLp3zscdPYH6Hl2GnEIDx_75MrOMpGNeYk)

I enjoyed writing this and it helped me get a better feel for n-body dynamics. Pull requests welcome :-)

## running the code

In this directory,
- run `cmake .`
- compile and run `driver.cc` with `make run`
- run the jupyter notebook to plot and animate the simulation trajectories

## requirements 

- C++: `cmake 3.11+` `gcc 10.2+` 
- python, jupyter, numpy, scipy, matplotlib, pandas
