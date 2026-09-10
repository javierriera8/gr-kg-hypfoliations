Compile the code with: 
$ make
To recompile, clear the previous compillation before, or the makefile will not see the changes: 
$ make clean; make

To run the code: 
exe/executable_name run/parameter_file_name.txt 

To run convergence tests (this feature may have to be updated), with  the number of runs (usually X=3): 
exe/executable_name run/parameter_file_name.txt -conv_test X 

Relevant parameters in parameter file: 
- file_prefix = 'folder_name/simulation_tag_': folder_name must exist, of the code will complain. 
... more to come ... 

ygraph or muninn (https://git.tpi.uni-jena.de/srenkhoff/muninn) can be used for visualization, although of course you can use anything you want. 
