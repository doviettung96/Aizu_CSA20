This assignment is about header files and make files. To do the assignment

Log in to envfs1 and do the following

1. copy or link the folders /home/csc09/makefile_demo and /home/csc09/makefile_assignment to your homefolder

2. the folder `makefile_demo` contains the demo makefile discussed in the class

3. edit the two files quad_serial.c and poisson_serial.c and

- break up each single file into separate files with one file per function

- break up each function file into header file (*.h) and implementation file (*.c)

- write header guards inside each header file.

- include header files in appropriate places using the #include macro

- finally write a makefile that can compile and link the executable as well as clean up the project, using make clean

4. Execute the output from each of the two original files and verify that your modifications provide the same results as the original executable.