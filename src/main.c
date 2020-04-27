#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
#include "consistency.h"
#include <time.h>

//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_circle_to_ellipse();
void dam_break();
void box();
void analyse_neighbours();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	 // script_circle_to_ellipse();
	//script_circle_to_ellipse();
	dam_break();
    //analyse_neighbours();
	return EXIT_SUCCESS;
}


void dam_break(){
	// Parameters of the problem
	double R = 0.1;
	double l = 0.057; // particle distribution on [-l,l] x [-l,l]
	double L = 1; // size of the domain: [-L,L] x [-L,L]
	double H = 1;
	double dt = 1.0e-4; // physical time step
	double T = 0.2; // duration of simulation
	bool gravity = 1; // 1 if we consider the gravity

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)


	// SPH parameters
	double N_x1 = 30;
	double N_y1 = 30;
	double N_x2 = 45;
	double N_y2 = 20;
	double kh = sqrt(21) * 2 * l / 25;

	//int n_per_dim = 60; // number of particles per dimension
	//double kh = sqrt(21) * 2 * l / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	
    int T_verlet = 4;
    double L_verlet = 1.1*(2*3.4*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,T_verlet);
    //void* verlet = NULL;
    
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	double CR = 1.0;
	double CF = 0.0;

	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	
	int n_p = N_x1 * N_y1 + N_x2 * N_y2;
	double h_x1 = (2 * l - 0.01) / (N_x1 - 1);
	double h_y1 = (2 * l - 0.01) / (N_y1 - 1);
	double h_x2 = (3 * l - 0.01) / (N_x2 - 1);
	double h_y2 = (0.5 * l ) / (N_y2 - 1);

	//int n_p = squared(n_per_dim);// +N_x * N_y; // total number of particles
	printf("nombre de point %d", n_p);
	//double h_x = (2 * l -0.01) / (n_per_dim - N_x2 - 1); // step between neighboring particles
	//double h_y= 2 * l / (n_per_dim - 1);
	//double m1 = rho_0 * h_x*h_y;
	double m1 = rho_0 * h_x1 * h_y1;
	double m2 = rho_0 * h_x2 * h_y2;
	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(n_p * sizeof(Residual*));

	for (int i = 0; i < N_y1; i++) {
		for (int j = 0; j < N_x1; j++) {
			int index = i * (N_x1)+j;
			xy* pos;
			xy* v;
			pos = xy_new(-2 * l + i * h_x1, j * h_y1);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m1, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	int index1 = N_x1 * N_y1;
	for (int i = 0; i < N_y2; i++) {
		for (int j = 0; j < N_x2; j++) {
			int index = index1 + i * (N_x2)+j;
			xy* pos;
			xy* v;
			pos = xy_new((j) * h_x2+0.005,-l+ i * h_y2);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m2, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	// Setup grid
	Grid *grid = Grid_new(-3*l, 4*l, -H, H, kh);
	// Setup BOUNDARY
	double lb = 0.420;
	double hb = 0.440;
	double Rp = 0.001; //particle radius
	Boundary* boundary = Boundary_new(-2*l-Rp,-Rp,3*l+Rp,-Rp,-l-Rp,hb-l+Rp,CR,CF);

	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
	// Simulation
	simulate_boundary(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation, boundary, &verlet);
	// Free memory
	Boundary_free(boundary);
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
/// with 2 liquids


// Evolution of a 2D circle with non-zero initial velocities (no surface tension force)
// Test case from "Simulating Free Surface Flows with SPH", Monaghan (1994)

