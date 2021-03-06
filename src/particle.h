#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "utils.h"

typedef struct Cell Cell;
typedef struct Grid Grid;
typedef struct Particle Particle;
typedef struct Particle_derivatives Particle_derivatives;
typedef struct Physical_parameters Physical_parameters;
typedef struct Verlet Verlet;


struct Cell {
	int i,j;
	List* neighboring_cells;
	List* particles;
	bool visited; // for improved algorithm
};

struct Grid {
	int nCellx;
	int nCelly;
	Cell*** cells;	// 2D array of cell pointers
	double h; // side length of a cell
	double left;	// x-coordinate of left side
	double right;	// x-coordinate of right side
	double top;		// y-coordinate of top side
	double bottom;	// y-coordinate of bottom side
};

struct Physical_parameters {
	double dynamic_viscosity;
	double rho_0;
	double gamma;
	double sound_speed;
	double sigma;
};

struct Particle {
	int index;
	double m;     // mass
	xy* pos;      // position
	xy* v;        // velocity
	double rho;   // density
	double P;     // pressure
	double Cs;    // color field
	xy* normal;   // normal vector
	double kappa; // curvature

	xy* XSPH_correction; // Correction on the velocity field when updating the position of the particles	
	bool on_free_surface; // boolean to know if particles is on the free surface (used for visualization)
	
	bool on_boundary; // boolean to detect if the particle is on the boundary and avoid the free surface force.
	
	Physical_parameters* param; // physical parameters associated to the particle

	Cell* cell;    // cell that the particle belongs to
	List* neighborhood; // list of neighbors
	List* potential_neighborhood; // list of potential neighbors (for Verlet)
};

struct Particle_derivatives {
	int index;
	double div_v;
	xy *grad_P;
	xy *lapl_v;
	xy *grad_Cs;
	xy *art_visc;
	double lapl_Cs;
};

struct Verlet {
	double L;
    double L_initial;
	int T;
    bool use_verlet;
};

// Grid
void Cell_free(Cell* cell); // Cell destructor

Grid* Grid_new(double x1, double x2, double y1, double y2, double kh); // Grid constructor
void Grid_free(Grid* grid); // Grid destructor

// Particle
Particle* Particle_new(int index, double m, xy* pos, xy* v, double rho_0, double mu, double c_0, double gamma, double sigma); // Particle constructor
void Particle_free(Particle* particle);
void free_particles(Particle** particles, int N);

double Particle_get_P(Particle *particle);
xy * Particle_get_v(Particle *particle);
xy * Particle_get_pos(Particle *particle);
double Particle_get_v_x(Particle *particle);
double Particle_get_v_y(Particle *particle);
double Particle_get_Cs(Particle *particle);
xy * Particle_get_normal(Particle *particle);

// Particle_derivatives
Particle_derivatives* Particle_derivatives_new(int index);
void Particle_derivatives_free(Particle_derivatives* particle_derivatives);
void Particle_derivatives_reset(Particle_derivatives *particle_derivatives);
void free_particles_derivatives(Particle_derivatives** particles_derivatives, int N);


// Update of the particles' locations in cells
Cell* localize_particle(Grid* grid, Particle * p);
void update_cells(Grid* grid, Particle** particles, int N);

// Update of the neighborhood
void add_neighbors_from_cell(Particle* p, Cell* cell, double r);
void add_neighbors_from_cells(Grid* grid, Particle* p);
void update_from_potential_neighbors(Particle** particles, int N, double r);
void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Verlet* verlet);

// Build random particles
Particle** build_particles(int N, double L);

void reset_grid(Grid* grid);
void reset_particles(Particle** particles, int N, int iter, Verlet* verlet);

// Verlet functions

Grid* Grid_new_verlet(double x1, double x2, double y1, double y2, double h, Verlet* v);

void Verlet_init(Verlet* v, double L, double L_initial, int T, bool use_verlet);

void update_verlet_cells(Grid* grid, Particle** particles, int N, Verlet* verlet);

Cell* localize_particle_verlet(Grid *grid, Particle *p, Verlet* verlet);

void create_potential_neighborhood(Grid* grid, Particle* p, Verlet* v);

void add_potential_neighbors_from_cell(Particle* p, Cell* cell , double r, double L);

#endif
