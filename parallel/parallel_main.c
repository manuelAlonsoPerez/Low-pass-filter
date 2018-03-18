#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

/*************************************************************
*  Matrix with the image values as a 2D-float array          * 
**************************************************************/
typedef struct{
    float** image_data;
    int m,n;
} image;


/*************************************************************
*  Allocates the struct image given the rows m and columns n * 
**************************************************************/
void allocate_image(image *u, int m, int n){
    int i;
    u->m = m;
    u->n = n;
    u->image_data = (float**)malloc(m*sizeof(float*));
    for (i = 0; i <m; i++){
        u-> image_data[i] = (float*)malloc(n*sizeof(float*));
    }

}

/*************************************************************
*  Deallocates the struct image                              * 
**************************************************************/
void deallocate_image(image *u){
    int i,n;

    n = sizeof (u->image_data);

    for ( i = 0; i < n; i++){
        free(u->image_data[i]);
    }
    free(u->image_data);
}

/***********************************************************************
*  Translates the values from the char array to the matrix             * 
***********************************************************************/
void convert_jpeg_to_image(const unsigned char* image_chars, image *u){
    int i,j;
    
    for(i = 0; i < u->m; i++){
        for(j = 0; j < u->n; j++){
            u->image_data[i][j]= (float)image_chars[i*u->n + j];
        }
    }
}

/***********************************************************************
*  Translates the values from  the matrix to the char array            * 
***********************************************************************/
void convert_image_to_jpeg( image* u, unsigned char *image_chars){
    
    int i,j;
    
    for(i = 0; i < u->m; i++){
        for(j = 0; j < u->n; j++){
            image_chars[i*u->n + j] =  (unsigned char) u->image_data[i][j];
        }
    }
    
}


/******************************************************************************************************************
*  Transfering of the rows in the borders of each partial image to the neighbouring processors subimages          * 
*******************************************************************************************************************/
void communicate2D_horizontal(image *u_local, int M_local, int N_local,int my_rank, int P){

    MPI_Status status;
    
    int lower_neigh = (my_rank>0) ? my_rank-1 : MPI_PROC_NULL;
    int upper_neigh = (my_rank<(P-1)) ? my_rank+1 : MPI_PROC_NULL;

    /* P0 case*/
    if(lower_neigh < 0 ){
        // send to upper neighbor, receive from upper neighor
        MPI_Sendrecv (&(u_local->image_data[M_local-2][0]),N_local,MPI_DOUBLE, upper_neigh,0,&(u_local->image_data[M_local-1][0]),
            N_local,MPI_DOUBLE, upper_neigh,MPI_ANY_TAG,MPI_COMM_WORLD, &status);
    }
    /* P-1 case*/
    else if(upper_neigh < 0 ){
        // send to lower neighbor, receive from lower neighor
        MPI_Sendrecv (&(u_local->image_data[1][0]),N_local,MPI_DOUBLE, lower_neigh,1,&(u_local->image_data[0][0]),
            N_local,MPI_DOUBLE, lower_neigh,MPI_ANY_TAG,MPI_COMM_WORLD, &status);
    }
    /* Intermediate P case*/
    else{

        // send to lower neighbor, receive from upper neighor
        MPI_Sendrecv (&(u_local->image_data[1][0]),N_local,MPI_DOUBLE, lower_neigh,2,&(u_local->image_data[M_local-1][0]),
            N_local,MPI_DOUBLE, upper_neigh,MPI_ANY_TAG,MPI_COMM_WORLD, &status);

        // send to upper neighbor, receive from lower neighor
        MPI_Sendrecv (&(u_local->image_data[M_local-2][0]),N_local,MPI_DOUBLE, upper_neigh,3, &(u_local->image_data[0][0]),
            N_local,MPI_DOUBLE, lower_neigh,MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    }    
}


int main (int argc, char *argv[]) {

    float kappa;
	int image_height,image_width, num_components,i, m, n, iters, iteration;
	int my_m, my_n, my_rank, num_procs, m_index, counter, counter_2;
	unsigned char *image_chars, *my_image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;
    int *sendcounts, *displs,*displs_2, *recvcounts;

    MPI_Status status;

	image u, u_bar,u_tmp;

    int mpi_rank = -1;
    int mpi_size = 0;

    /* Inicializing the MPI environment */
    MPI_Init(&argc, &argv);
    
    /* Reading the number on processors and  assigning each one his number */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    /*------------------------------------------------------------------------------------- 
    -  Reads from command line kappa, the number of iterations, the input and outpufiles. -                
    -  If no paramameters are passed a fix set is defined.                                -
    -------------------------------------------------------------------------------------*/
    if(argc >= 5){
        kappa = atof(argv[1]);
        iters = atoi(argv[2]);
        input_jpeg_filename = argv[3];
        output_jpeg_filename = argv[4]; 
    } else {
        kappa = 0.1;
        iters = 50;
        input_jpeg_filename = "mona_lisa_noisy.jpeg";
        output_jpeg_filename = "mona_lisa_parallel_k=0.1_it=50.jpeg";
    }

    
    /* The first processor reads the image file and assigns the corresponding values to image_chars, 
    image_height, image_width, num_components  */

	if (my_rank==0) {

		import_JPEG_file (input_jpeg_filename, &image_chars,&image_height,&image_width,&num_components);

        printf("Kappa is: %f\n", kappa);
        printf("Total iterations: %d\n", iters);
        printf("The input file is : %s\n", input_jpeg_filename);
        printf("The output file is : %s\n", output_jpeg_filename);
	}

    /* Sends to each processor image_height and image_weight */
    MPI_Bcast (&image_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&image_width, 1, MPI_INT, 0, MPI_COMM_WORLD);



    /*----------------------------------------------------------------------------------------------------- *
    - Dividing the work, each process will have a portion of the array image_chars, and  will construct the -
    - corresponding image struct to do the computations. Each processor will have and extra row two make    -
    - possible the computing in the border areas, two extra rows in the case of the intermediate processes. -                                            
    /*---------------------------------------------------------------------------------------------------- */


    /* The array sedcounts  records how many chars each process will recieve, and the array dspls 
    indicates the index of image_chars to start sending elements*/

    sendcounts = malloc(num_procs*sizeof(int));
    displs = malloc(num_procs*sizeof(int));
    counter = 0;
    
    my_n = image_width;

    for(i = 0; i < num_procs; i++) {
        m_index = ((i + 1)*image_height)/num_procs - (( i * image_height ) / num_procs);
        
        if(i == 0 || i == num_procs-1){sendcounts[i] = (m_index +1) * my_n ;}
        else{sendcounts[i] = (m_index + 2) * my_n ;}
        
        displs[i] = counter;
        
        if(i == 0){counter += (m_index-1) * my_n;}
        else{counter += (m_index) * my_n;}
    }

    /* Allocating the local char array my_image_chars and  the local images u_, u_bar corresponding to the process */
    
    my_image_chars = malloc(sendcounts[my_rank]*image_width*sizeof(unsigned char));
    my_m = sendcounts[my_rank]/my_n;
   
    allocate_image (&u, my_m, my_n);
    allocate_image (&u_bar, my_m, my_n);
    
    //printf("\nI am process number %d, My rows are: %d, My columns  are: %d \n",my_rank, sendcounts[my_rank]/my_n, my_n );

    /* each process asks process 0 for a partitioned region of image_chars and copy the values into u */
    
    MPI_Scatterv(&image_chars[0], sendcounts, displs, MPI_UNSIGNED_CHAR,&my_image_chars[0], sendcounts[my_rank] * my_n, MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
    convert_jpeg_to_image(my_image_chars, &u);

    /*------------------------------------------------------------------------------------------------------------ *
    - Doing the computing, each processor performs the denoising algorithm on his partial matrix for the specified -
    - number of itereations.On each iteration the neigbouring processes exchange the first and last rows computed  -
    /*----------------------------------------------------------------------------------------------------------- */

    iteration = 0;

    while(iteration  < iters ){
        
        
        for(m = 1; m < my_m-1; m++){
            for(n = 1; n < my_n-1; n++){
                u_bar.image_data[m][n] = u.image_data[m][n] + (kappa*(u.image_data[m-1][n] +
                                          u.image_data[m][n-1]- (4*u.image_data[m][n]) +                                                                               
                                          u.image_data[m][n+1]+ u.image_data[m+1][n])) ;
            }
        }

        if (my_rank == 0 || iteration == 0){
            for(n = 0; n < my_n; n++){
                u_bar.image_data[0][n] = u.image_data[0][n];
            }

        }
        
        if (my_rank == num_procs-1 || iteration == 0){
            for(n = 0; n < my_n; n++){   
                u_bar.image_data[my_m-1][n] = u.image_data[my_m-1][n];
            }

        }
        
        for(m = 0; m < my_m; m++){
            u_bar.image_data[m][0] = u.image_data[m][0];
            u_bar.image_data[m][my_n-1] = u.image_data[m][my_n-1];
        }
        
        communicate2D_horizontal(&u_bar, u_bar.m, u_bar.n,my_rank, num_procs);

        iteration++;
        u_tmp = u;
        u = u_bar;
        u_bar = u_tmp;
    }
    
    
    /*---------------------------------------------------------------------------------------------------------------- *
    - Each process sends its resulting content of u_bar, transformned in a char array again, to process 0.             -
    - Process 0 receives from each process the incoming values and copy them into the designated region of image_chars -                                     -
    /*---------------------------------------------------------------------------------------------------------------  */
    
    convert_image_to_jpeg(&u_bar, my_image_chars);
    deallocate_image (&u);
    deallocate_image (&u_bar);
    
    /* The array recvcounts  records how many chars each process will recieve, and the array dspls_2 
    indicates the index of image_chars to start sending elements*/
    recvcounts = (int *)malloc(num_procs*sizeof(int));
    displs_2 = malloc(num_procs*sizeof(int));
    counter_2 = 0;

    for(i = 0; i < num_procs; i++) {
        if(i == 0){recvcounts[i] = sendcounts[i] - my_n ;}
        else if(i == num_procs-1){recvcounts[i] = sendcounts[i] - my_n ;}
        else{ recvcounts[i] = sendcounts[i] - 2*my_n ;}
        displs_2[i] = counter_2;
        counter_2 += recvcounts[i] ;
    }

    //printf("\nI am process number %d, I send chars from index: %d to index: %d\n",my_rank, displs_2[my_rank], displs_2[my_rank]+recvcounts[my_rank]) ;


    MPI_Gatherv(&my_image_chars[0], recvcounts[my_rank],  MPI_UNSIGNED_CHAR, &image_chars[-1], recvcounts,displs, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    /* The processor 0 converts the char array image_chars to and image file and deallocates the whole image */
    if (my_rank==0) {
        export_JPEG_file (output_jpeg_filename, image_chars,image_height, image_width,num_components, 75);
    }

    MPI_Finalize();



  return 0;
}
