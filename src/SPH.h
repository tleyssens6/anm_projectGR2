#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"

typedef struct Boundary Boundary;
typedef struct Setup Setup;
typedef struct Residual Residual;
typedef enum Free_surface_detection Free_surface_detection;


enum Free_surface_detection {CSF,DIVERGENCE,NONE};

struct Setup{
	int itermax;
	double timestep;
	double kh;
	Verlet* verlet;
	Kernel kernel;
	Free_surface_detection free_surface_detection;
	double interface_threshold;
	double XSPH_epsilon;
	bool gravity;
};

struct Residual {
	double mass_eq;
	double momentum_x_eq;
	double momentum_y_eq;
};

struct Boundary{
  double xleft1;
  double xleft2;
  double xright;
  double ybottom1;
  double ybottom2;
  double ytop;
	double CR;
  double CF;
};

Setup* Setup_new(int iter, double timestep,double kh, Verlet* verlet, Kernel kernel,Free_surface_detection free_surface_detection,double interface_threshold, double XSPH_epsilon, bool gravity);
void Setup_free(Setup* setup);

double max_velocity(Particle** p, int n_p);

Residual* Residual_new();
void free_Residuals(Residual** residuals,int N);

void simulate_SPH(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup, Animation* animation, Boundary* boundary);
void update_positions(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);

void compute_Cs(Particle *particle, Kernel kernel, double kh);
void assemble_residual_NS(Particle* particle, Particle_derivatives* Particle_derivatives, Residual* residual, Setup* setup);
void time_integrate(Particle* particle, Residual* residual, double delta_t);

void compute_XSPH_correction(Particle *particle, Kernel kernel, double kh,double epsilon);

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives);

Boundary* Boundary_new(double xleft1, double xleft2, double xright, double ybottom1, double ybottom2, double ytop, double CR, double CF);
void Boundary_free(Boundary* boundary);
void reflective_boundary(Particle** p, int n_p, Boundary* boundary,double Rp);

xy* compute_surfaceTension(Particle* p, Particle_derivatives* d);

#endif
