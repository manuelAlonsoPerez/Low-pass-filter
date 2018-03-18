#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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

/*************************************************************************************************
*  Performs the denoising algorithm for the matrix given, accepts parameter kappa, iters         * 
**************************************************************************************************/
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters){

    int iteration, m,n;
    
    image * u_tmp;
    iteration = 0;
    
    while(iteration  < iters ){
        /* Performing the computations in the inner points of the matrix */
        for(m = 1; m < u->m-1; m++){
            for(n = 1; n  < u->n-1; n++){
                u_bar->image_data[m][n] = u->image_data[m][n] + (kappa*(u->image_data[m-1][n] +
                                          u->image_data[m][n-1]- (4*u->image_data[m][n]) +                                                                               u->image_data[m][n+1]+ u->image_data[m+1][n])) ;
            }
        }
        /* Copying the first and last row*/
        for(n = 0; n < u->n; n++){
            u_bar->image_data[0][n] = u->image_data[0][n];
            u_bar->image_data[u->m-1][n] = u->image_data[u->m-1][n];
        }
        /* Copying the first and last Column*/
        for(m = 0; m < u->m; m++){
            u_bar->image_data[m][0] = u->image_data[m][0];
            u_bar->image_data[m][u->n-1] = u->image_data[m][u->n-1];
        }
        /* Advancing to next iteration and assigning u the values of u_bar for the next iteration*/
        iteration++;
        u_tmp = u;
        u = u_bar;
        u_bar = u_tmp;
        
    }
}


int main (int argc, char *argv[]) {
    
    float kappa;
	int image_height,image_width,num_components,i, iters;
	unsigned char *image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;

	image u,u_bar;
    
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
        output_jpeg_filename = "mona_lisa_serial_k=0.1_it=50.jpeg";
    }
    
    printf("kappa is: %f\n", kappa);
    printf("iters is: %d\n", iters);
    printf("the input file is : %s\n", input_jpeg_filename);
    printf("the output file is : %s\n", output_jpeg_filename);

    /* Import the image specified, assigning the corresponding values to image_chars, image_height, image_width, num_components  */
    import_JPEG_file (input_jpeg_filename, &image_chars,&image_height,&image_width,&num_components);
    printf("Image Name: %s ,Image Width:%d , Image Height:%d , Components:%d\n","mona_lisa_noisy.jpg", image_width, image_height,num_components);
     
    allocate_image(&u ,image_height,image_width);       /* Allocating the two images ussed for the computation*/
    allocate_image(&u_bar ,image_height,image_width);
    
    
    convert_jpeg_to_image(image_chars,&u);  /* Translating the char array values to the image matrix values*/
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);  /* Performing the denoising algorithm for the values of kappa and iterations given  */
    convert_image_to_jpeg(&u_bar, image_chars);  /* Converting back the matrix values to char array  */
    export_JPEG_file (output_jpeg_filename, image_chars,image_height, image_width,num_components, 75); /* Creating the image from the array */
    
    deallocate_image(&u);      /* Deallocating the two images ussed for the computation*/
    deallocate_image(&u_bar);
    
	return 0;    
    
}
