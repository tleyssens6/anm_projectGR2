#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
#include <time.h>

void drop_no_pool();
void dam_break();
void drop_with_pool();



int main() {
    //dam_break();
	//drop_no_pool();
    drop_with_pool();
	return EXIT_SUCCESS;
}

void dam_break() {
    
    /* Dimensional Parameters of the problem */
    double drop_height = 0.0;           // height of the drop : should be between 0 and 1
    double L = 0.5;                     // size of the domain: [-L,L] x [-L,L]
    double dt = 5.e-5;                  // physical time step
    double T = 5.;                      // duration of simulation
    int n_iter = (int)(T / dt);         // number of iterations to perform
    
    
    /* Physical parameters */
    bool gravity = true;
    double rho_0 = 1000.0;              // flow
    double mu_0 = 1.0016e-3;            // dynamic viscosity of water at 20°C (in N.s/m^2)
    double gamma = 7.0;                 // typical value for liquid (dimensionless)
    double c_0 = 1.0;                   // sound speed in water at 20°C (in m/s)
    double sigma = 72.86e-3;            // surface tension of water-air interface at 20°C (in N/m)
    double XSPH_epsilon = 0.5;          // parameter used in the XSPH correction
    Free_surface_detection surface_detection = DIVERGENCE; // Surface detection using the divergence of the position
    double interface_threshold = 1.5;   // threshold to determine if particle is on free surface
    
    
    /* Initialization of particles */
    int N_x1 = 61;                      // number of particles in x direction
    int N_y1 = 61;                      // number of particles in y direction
    int n_p = N_x1 * N_y1;
    printf("number of particles : %d\n", n_p);
    double dx_initial = (0.8*L) / (N_x1); // initial distance between 2 particles
    printf("dx_initial = %f\n", dx_initial);
    double kh = 3.1*dx_initial;           // smoothing length
    printf("kh = %f\n",kh);
    Kernel kernel = Cubic;              // kernel choice

    double m1 = rho_0 * dx_initial * dx_initial;
    
    Particle** particles = (Particle**)malloc((n_p) * sizeof(Particle*));
    Particle_derivatives** particles_derivatives = malloc((n_p) * sizeof(Particle_derivatives*));
    Residual** residuals = malloc((n_p) * sizeof(Residual*));
    
    for (int i = 0; i < N_x1; i++) {
        for (int j = 0; j < N_y1; j++) {
            int index = i*N_y1 + j;
            xy* pos;
            xy* v;
            pos = xy_new(-L + i * dx_initial+0.0005*L,(-1+drop_height)*L + j * dx_initial+0.0005*L);
            v = xy_new(0.0, 0.0);
            particles[index] = Particle_new(index, m1, pos, v, rho_0, mu_0, c_0, gamma, sigma);
            particles_derivatives[index] = Particle_derivatives_new(index);
            residuals[index] = Residual_new();
        }
    }
    
    int N_flow = N_x1*N_y1;
    
    
    /* Setup grid using verlet or not */
    int T_verlet = 4;
    double L_verlet = 1.1*(2*3.*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    double L_verlet_initial = L_verlet;
    printf("L_Verlet = %f \n",L_verlet);
    bool use_verlet = true;
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,L_verlet_initial, T_verlet,use_verlet);
    Grid* grid;
    if (use_verlet) {
        grid = Grid_new_verlet(-1.4*L, 1.4*L,-1.4*L, 1.4*L , kh, &verlet);
    }
    else {
        grid = Grid_new(-2.*L, 2.*L,-2.*L, 2.*L , kh);
    }
    
    
    /* Setup BOUNDARY */
    double CR = 0.7;                                // reflection coefficient
    double CF = 0.3;                                // friction coefficient
    Boundary* boundary = Boundary_new(-L, -0.2*L, L+0.001, (-1+drop_height)*L, -L, L, CR, CF);

    
    /* Setup setup */
    Setup* setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
    
    
    /* Setup animation */
    double T_anim = 10;                             // duration of animation
    double dt_anim = T_anim / n_iter;               // time step of animation
    int N_pool = 0;
    Animation* animation = Animation_new(N_pool, N_flow, dt_anim,grid, 1);
    
    
    /* Simulation */
    simulate_SPH(grid, particles, particles_derivatives, residuals, n_p, setup, animation, boundary);
    
    
    /* memroy freeing */
    Boundary_free(boundary);
    free_particles(particles, n_p);
    free_particles_derivatives(particles_derivatives, n_p);
    free_Residuals(residuals, n_p);
    Grid_free(grid);
    Setup_free(setup);
    Animation_free(animation);
}



void drop_no_pool() {
    
    /* Dimensional Parameters of the problem */
    double drop_height = 0.6;           // height of the drop : should be between 0 and 1
    double L = 0.5;                     // size of the domain: [-L,L] x [-L,L]
	double dt = 5.e-5;                  // physical time step
    double T = 5.;                      // duration of simulation
    int n_iter = (int)(T / dt);         // number of iterations to perform
    
    
	/* Physical parameters */
    bool gravity = true;
	double rho_0 = 1000.0;              // flow
	double mu_0 = 1.0016e-3;            // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0;                 // typical value for liquid (dimensionless)
	double c_0 = 1.0;                   // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3;            // surface tension of water-air interface at 20°C (in N/m)
    double XSPH_epsilon = 0.5;          // parameter used in the XSPH correction
    Free_surface_detection surface_detection = DIVERGENCE; // Surface detection using the divergence of the position
    double interface_threshold = 1.5;   // threshold to determine if particle is on free surface
    
    
    /* Initialization of particles */
    int N_x1 = 81;                      // number of particles in x direction
    int N_y1 = 41;                      // number of particles in y direction
    int n_p = N_x1 * N_y1;
    printf("number of particles : %d\n", n_p);
    double dx_initial = (0.8*L) / (N_x1); // initial distance between 2 particles
    printf("dx_initial = %f\n", dx_initial);
    double kh = 3.1*dx_initial;           // smoothing length
    printf("kh = %f\n",kh);
    Kernel kernel = Cubic;              // kernel choice

	double m1 = rho_0 * dx_initial * dx_initial;
    
	Particle** particles = (Particle**)malloc((n_p) * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc((n_p) * sizeof(Particle_derivatives*));
	Residual** residuals = malloc((n_p) * sizeof(Residual*));
    
	for (int i = 0; i < N_x1; i++) {
		for (int j = 0; j < N_y1; j++) {
			int index = i*N_y1 + j;
			xy* pos;
			xy* v;
            pos = xy_new(-L + i * dx_initial+0.0005*L,(-1+drop_height)*L + j * dx_initial+0.0005*L);
			v = xy_new(0.0, 0.0);
			particles[index] = Particle_new(index, m1, pos, v, rho_0, mu_0, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
    
    int N_flow = N_x1*N_y1;
    
    
	/* Setup grid using verlet or not */
    int T_verlet = 4;
    double L_verlet = 1.1*(2*3.*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    double L_verlet_initial = L_verlet;
    printf("L_Verlet = %f \n",L_verlet);
    bool use_verlet = true;
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,L_verlet_initial, T_verlet,use_verlet);
    Grid* grid;
    if (use_verlet) {
        grid = Grid_new_verlet(-1.4*L, 1.4*L,-1.4*L, 1.4*L , kh, &verlet);
    }
    else {
        grid = Grid_new(-2.*L, 2.*L,-2.*L, 2.*L , kh);
    }
    
    
	/* Setup BOUNDARY */
    double CR = 0.7;                                // reflection coefficient
    double CF = 0.3;                                // friction coefficient
    Boundary* boundary = Boundary_new(-L, -0.2*L, L, (-1+drop_height)*L, -L, L, CR, CF);

    
	/* Setup setup */
	Setup* setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
    
    
	/* Setup animation */
    double T_anim = 10;                             // duration of animation
    double dt_anim = T_anim / n_iter;               // time step of animation
    int N_pool = 0;
    Animation* animation = Animation_new(N_pool, N_flow, dt_anim,grid, 1);
    
    
	/* Simulation */
	simulate_SPH(grid, particles, particles_derivatives, residuals, n_p, setup, animation, boundary);
    
	
    /* memroy freeing */
	Boundary_free(boundary);
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

void drop_with_pool() {
    /* Dimensional Parameters of the problem */
    double drop_height = 0.4;           // height of the drop : should be between 0 and 1
    double L = 0.5;                     // size of the domain: [-L,L] x [-L,L]
    double dt = 5.e-5;                  // physical time step
    double T = 5.;                      // duration of simulation
    int n_iter = (int)(T / dt);         // number of iterations to perform
    
    
    /* Physical parameters */
    bool gravity = true;
    double rho_0 = 1000.0;              // flow
    double rho_1 = 800.;                // pool
    double mu_0 = 1.0016e-3;            // dynamic viscosity of water at 20°C (in N.s/m^2)
    double mu_1 = 0.086;                // pool dynamic viscosity
    double gamma = 7.0;                 // typical value for liquid (dimensionless)
    double c_0 = 1.0;                   // sound speed in water at 20°C (in m/s)
    double sigma = 72.86e-3;            // surface tension of water-air interface at 20°C (in N/m)
    double XSPH_epsilon = 0.5;          // parameter used in the XSPH correction
    Free_surface_detection surface_detection = DIVERGENCE; // Surface detection using the divergence of the position
    double interface_threshold = 1.5;   // threshold to determine if particle is on free surface
    
    
    /* Initialization of particles */
    int N_x1 = 51;                      // number of particles in x direction
    int N_y1 = 41;                      // number of particles in y direction
    int N_x2 = ceil(N_x1*3/2)+1;        // number of particles in x direction for the pool
    int N_y2 = ceil(N_x2/6);            // number of particles in y direction for the pool
    int n_p = N_x1*N_y1 + N_x2*N_y2;
    printf("number of particles : %d\n", n_p);
    double dx_initial = (0.8*L) / (N_x1); // initial distance between 2 particles
    printf("dx_initial = %f\n", dx_initial);
    double kh = 3.1*dx_initial;           // smoothing length
    printf("kh = %f\n",kh);
    Kernel kernel = Cubic;              // kernel choice

    double m1 = rho_0 * dx_initial * dx_initial;
    double m2 = rho_1 * dx_initial * dx_initial;
    
    Particle** particles = (Particle**)malloc((n_p) * sizeof(Particle*));
    Particle_derivatives** particles_derivatives = malloc((n_p) * sizeof(Particle_derivatives*));
    Residual** residuals = malloc((n_p) * sizeof(Residual*));
    
    for (int i = 0; i < N_x1; i++) {
        for (int j = 0; j < N_y1; j++) {
            int index = i*N_y1 + j;
            xy* pos;
            xy* v;
            pos = xy_new(-L + i * dx_initial+0.0005*L,(-1+drop_height)*L + j * dx_initial+0.0005*L);
            v = xy_new(0.0, 0.0);
            particles[index] = Particle_new(index, m1, pos, v, rho_0, mu_0, c_0, gamma, sigma);
            particles_derivatives[index] = Particle_derivatives_new(index);
            residuals[index] = Residual_new();
        }
    }
    int index0 = N_x1*N_y1;
    for (int i = 0; i < N_x2; i++) {
        for (int j = 0; j < N_y2; j++) {
            int index = index0 + i * (N_y2)+j;
            xy* pos;
            xy* v;
            pos = xy_new(-0.2*L + i*dx_initial +0.0005*L, -L + j * dx_initial);
            v = xy_new(0.0, 0.0); // initial velocity = 0
            particles[index] = Particle_new(index, m2, pos, v, rho_1, mu_1, c_0, gamma, sigma);
            particles_derivatives[index] = Particle_derivatives_new(index);
            residuals[index] = Residual_new();
        }
    }
    int N_flow = N_x1*N_y1;
    int N_pool = N_x2*N_y2;
    
    /* Setup grid using verlet or not */
    int T_verlet = 4;
    double L_verlet = 1.1*(2*3.*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    double L_verlet_initial = L_verlet;
    printf("L_Verlet = %f \n",L_verlet);
    bool use_verlet = true;
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,L_verlet_initial, T_verlet,use_verlet);
    Grid* grid;
    if (use_verlet) {
        grid = Grid_new_verlet(-1.4*L, 1.4*L,-1.4*L, 1.4*L , kh, &verlet);
    }
    else {
        grid = Grid_new(-2.*L, 2.*L,-2.*L, 2.*L , kh);
    }
    
    
    /* Setup BOUNDARY */
    double CR = 0.7;                                // reflection coefficient
    double CF = 0.3;                                // friction coefficient
    Boundary* boundary = Boundary_new(-L, -0.2*L, L, (-1+drop_height)*L, -L, L, CR, CF);

    
    /* Setup setup */
    Setup* setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
    
    
    /* Setup animation */
    double T_anim = 10;                             // duration of animation
    double dt_anim = T_anim / n_iter;               // time step of animation
    Animation* animation = Animation_new(N_pool, N_flow, dt_anim,grid, 1);
    
    
    /* Simulation */
    simulate_SPH(grid, particles, particles_derivatives, residuals, n_p, setup, animation, boundary);
    
    
    /* memroy freeing */
    Boundary_free(boundary);
    free_particles(particles, n_p);
    free_particles_derivatives(particles_derivatives, n_p);
    free_Residuals(residuals, n_p);
    Grid_free(grid);
    Setup_free(setup);
    Animation_free(animation);
}

