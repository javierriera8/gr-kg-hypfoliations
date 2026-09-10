Compile the code with: 
$ make
To recompile, clear the previous compillation before, or the makefile will not see the changes: 
$ make clean; make

To run the code: 
exe/executable_name run/parameter_file_name.txt 

To run convergence tests (this feature may have to be updated), with  the number of runs (usually X=3): 
exe/executable_name run/parameter_file_name.txt -conv_test X 

Relevant parameters in parameter file: 

* File creation:
  - file_prefix = 'folder_name/simulation_tag_': folder_name must exist, of the code will complain. 
* Physical parameters:
  - parameters%physics%gr: used to control the use (=0) or not (=1) of the Cowling approximation
  - parameters%physics%masskgf: massive term at the Klein-Gordon field
  - parameters%physics%a, sigma, center: control of the initial data for the real scalar field
  - parameters%physics%aim, sigmaim, centerim: control of the initial data for the imaginary scalar field
* Numerical parameters:
  - parameters%timer%end_time: ending time of the simulation, for example=20
  - parameters%timer%dt: width of the time cell, for example=0.001d0 
  - parameters%grid%ncells: number of cells for the spatial grid, for example = 200
  - parameters%output%every_t_step: number of time-steps before the input is written, for example = 200 

To visualize general simulations, ygraph or muninn (https://git.tpi.uni-jena.de/srenkhoff/muninn) is recommendend. To visualize convergence tests, the notebooks included are recommended. Nevertheless, of course you can use anything you want.
