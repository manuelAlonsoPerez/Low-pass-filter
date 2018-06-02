## PROJECT DESCRIPTION

Noise removing of an image preforming a linnear lowpass filter, the filter 
operations are performed parallel in the number of cores selected when 
compiling.


## COMPILING INSTRUCTIONS       

In linux terminal: 

1.- Execute Makefile in main folder to compile the serial and parallel programs 
and to compile all the functions and make the lib-simple library	


`$ make all`


2.-  To execute serial_main in linux terminal
	`$ ./serial_main  `//automatic execution with kappa=0.1 and 50 iters
	`$ ./serial_main kappa iters infile outfile `// Passing parameters


3.-  To execute serial_main in abel
	`$ sbatch serial_main.scp ` //automatic execution with kappa=0.1 and 50 iters
   `$ sbatch serial_main.scp kappa iters infile outfile` // Passing parameters


4.-  To execute parallel_main in linux terminal
	`$ mpirun -n nr_processors ./parallel_main ` //automatic execution  kappa=0.1 and 50 iters
	`$ mpirun -n nr_processors ./parallel_main kappa iters infile outfile` //Passing parameters


5.-  To execute parallel_main in abel
`$ sbatch parallel_main.scp`  //automatic execution with kappa=0.1 and 50 iters
`$ sbatch parallel_main.scp kappa iters infile outfile` // Passing parameters

## PACKAGE CONTENT        

**Makefile: **
Compiles the serial and parallel programs and to compile all the functions and make  the lib-simple library, remove the objects and reports from execution on request and compress all for delivery.

**serial:**
Makefile that compiles the serial version and the simple-jpeg library. serial_main.c the code of the serial version of the assignment serial_main.scp  the script to run the serial version on abel mona_lisa_noisy.jpeg the image use for denoising. 

***.jpeg  ** the result of running the program on abel for different kappa values and iterations

**parallel:**
Makefile that compiles the parallel version and the simple-jpeg library.
parallel_main.c the code of the parallel version of the assignment
parallel_main.scp  the script to run the parallel version on abel
mona_lisa_noisy.jpeg the image use for denoising.
***.jpeg  **the result of running the program on abel for different kappa values and iterations


**simple-jpeg:**
A set of functions to work on jpeg images, we will use import_JPEG_file and export_JPEG_file which are in a build library.


**Readme.txt:**
this

## PECULIARITIES          

In the assignment text it is suggested the use of a whole image  and a sub image for each process. I have used a whole char array and a partial array for each process. Building this partial images from the corresponding portion of image_chars. The image is divided horizontally trying to assign each processor similar amount of information.
The first and last process take one extra row to allow computation in the border row. This last computed row is send to the neighbor processes. The processes in the middle take one extra upper and downer row for the same reason, and send the two last rows  computed to the neighbor processes.
After each processor has finished the computation converts the image to a char array and sends it to the first process where the array is gathered and the whole denoised image built back.
 

## PROGRAM   STATUS       


After testing on abel with one core and 8 processors and the scripts provided, both  serial and parallel version works with and without passing arguments, the level of  detail in the images seems to concord with the values given.
