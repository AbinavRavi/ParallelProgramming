#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <pthread.h>
#include <math.h>

#include "mandelbrot_set.h"
#define CHUNK_SIZE 10

int chunksize;
int chunkcount;
int current_chunk;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct
{
	int x_resolution;
	int y_resolution;
	int max_iter;
	double view_x0;
	double view_x1;
	double view_y0;
	double view_y1;
	double x_stepsize;
	double y_stepsize;
	int palette_shift;
	void *image;
	int num_threads;
}pthread_args;

void *compute_mandelbrot(void *args)
{
		pthread_args *args1 = args;

		pthread_mutex_lock(&mutex);

		int task_id = current_chunk;
		if (current_chunk == chunkcount){
			pthread_mutex_unlock(&mutex);
			return NULL;
		}
		chunkcount += chunksize;
		pthread_mutex_unlock(&mutex);

		int row_start = chunksize;
		int row_end = row_start+ chunksize;
		if (task_id == chunkcount - 1){
			row_end = args1->y_resolution;
		}

		double y;
		double x;

	unsigned char (*img)[args1->x_resolution][3]=(unsigned char(*)[args1->x_resolution][3])args1->image;

	complex double Z;
	complex double C;

	int k;

	for (int i = row_start; i < row_end; i++)
		{
			for (int j = 0; j < args1->x_resolution; j++)
			{
				y = args1->view_y1 - i * args1->y_stepsize;
				x = args1->view_x0 + j * args1->x_stepsize;

				Z = 0 + 0 * I;
				C = x + y * I;

				k = 0;

				do
				{
					Z = Z * Z + C;
					k++;
				} while (cabs(Z) < 2 && k < (args1->max_iter));

				if (k == args1->max_iter)
				{
					memcpy(img[i][j], "\0\0\0", 3);
				}
				else
				{
					int index = (k + args1->palette_shift)
					% (sizeof(colors) / sizeof(colors[0]));
					memcpy(img[i][j], colors[index], 3);
				}
			}
		}


}

void mandelbrot_draw(int x_resolution, int y_resolution, int max_iter,double view_x0, double view_x1, double view_y0, double view_y1,double x_stepsize, double y_stepsize,int palette_shift, unsigned char (*image)[x_resolution][3],int num_threads) {

						pthread_t *threads = (pthread_t*) malloc (num_threads*sizeof(pthread_t));
						pthread_args* args = (pthread_args*)malloc(num_threads*sizeof(pthread_args));


							for (int z = 0; z<num_threads;++z)
							{
								args[z].num_threads = num_threads;
								args[z].x_resolution = x_resolution;
								args[z].y_resolution = y_resolution;
								args[z].max_iter = max_iter;
								args[z].view_x0 = view_x0;
								args[z].view_x1 = view_x1;
								args[z].view_y0 = view_y0;
								args[z].view_y1 = view_y1;
								args[z].image = (void*) image;
								args[z].palette_shift = palette_shift;
								args[z].x_stepsize = x_stepsize;
								args[z].y_stepsize = y_stepsize;
								pthread_create(&threads[z],NULL, compute_mandelbrot ,&args[z]);
							}

							chunksize = floor(y_resolution /num_threads);
							chunksize = chunksize > CHUNK_SIZE? CHUNK_SIZE :chunksize;
							chunkcount = ceil(y_resolution/chunksize);

							current_chunk = 0;

							for (int i = 0;i<num_threads;++i)
							{
							pthread_join(threads[i],NULL);
							}

							free(threads);
							free(args);



}
