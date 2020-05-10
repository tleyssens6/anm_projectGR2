#ifndef PRINT_PARTICULES_H
#define PRINT_PARTICULES_H

#include "BOV.h"
#include <math.h>
#include "particle.h"
#include "utils.h"


typedef struct Animation Animation;

struct Animation {
	bov_window_t* window;
	bov_points_t* particles_pool;
    bov_points_t* particles_flow;
	double timeout;
	int N_pool;
    int N_flow;
	bov_points_t* grid;
};

Animation* Animation_new(int N_pool, int N_flow, double timeout,Grid* grid,double scale);
void Animation_free(Animation* animation);

void fillData(GLfloat(*data_pool)[8], GLfloat(*data_flow)[8], Particle** particles, int N_pool, int N_flow);
bov_points_t * load_Grid(Grid* grid,double scale);
void display_particles_boundary(Particle** particles, Animation* animation,bool end, int iter, double bounds[6]);

// thomas functions

void display_neighbours(bov_window_t* window, Animation* animation, Particle** particles, int N);
void fillData_pressureGrad(GLfloat(*data_pool)[8], GLfloat(*data_flow)[8], Particle** particles, int N_pool, int N_flow);



#endif
