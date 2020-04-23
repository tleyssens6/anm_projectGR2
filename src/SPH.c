#include "SPH.h"
#include <time.h>

Setup* Setup_new(int iter, double timestep,double kh,Verlet* verlet,Kernel kernel, Free_surface_detection free_surface_detection, double interface_threshold,double XSPH_epsilon, bool gravity) {
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->kh = kh;
	setup->verlet = verlet;
	setup->kernel = kernel;
	setup->free_surface_detection = free_surface_detection;
	setup->interface_threshold = interface_threshold;
	setup->XSPH_epsilon = XSPH_epsilon;
	setup->gravity = gravity;
	return setup;
}

void Setup_free(Setup* setup) {
	if(setup->verlet != NULL)
		free(setup->verlet);
	free(setup);
}

Residual* Residual_new(){
	Residual* residual = (Residual*)malloc(sizeof(Residual));
	residual->mass_eq = 0;
	residual->momentum_x_eq = 0;
	residual->momentum_y_eq = 0;
	return residual;
}

void free_Residuals(Residual** residuals, int N) {
	for (int i = 0;i < N;i++)
		free(residuals[i]);
	free(residuals);
}
double max_velocity(Particle** p, int n_p){
	double max = 0;
	for(int i = 0 ; i < n_p; i++){
		double v = fabs(p[i]->v->x);
		if(max < v){
			max = v;
		}
	}
	return max;
}

void simulate_boundary(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation, Boundary* boundary, Verlet* verlet){
    double current_time = 0.0;
    double Rp = 0.001;
    int ii = 5;
    double bounds[4] = {boundary->xleft1,boundary->xright, boundary->ybottom2, boundary->ytop};
    printf("%d\n", setup->itermax);
    for (int iter = 0; iter < setup->itermax; iter++) {
        printf("----------------------------------------------------- \n");
        printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
        //clock_t t0 = clock();
        update_verlet_cells(grid, particles, n_p, verlet);
        update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
        //clock_t t1 = clock();
        //printf("time for neighborhood update without verlet : %f \n",(double)(t1-t0));
        if (animation != NULL){
            display_particles_boundary(particles, animation, false,iter,bounds);
        }
        update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);

        printf("velocity_max = %f\n", max_velocity(particles,n_p));
        reflective_boundary(particles,n_p,boundary,Rp);
        if (iter%ii == 0){
            // density_correction_MLS(particles, n_p, setup->kh, setup->kernel);
        }
        get_M0(particles,n_p,setup->kh,setup->kernel);
        get_M1(particles,n_p,setup->kh,setup->kernel);
        current_time += setup->timestep;
    }
    
    update_verlet_cells(grid, particles, n_p, verlet);
    update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
    if (animation != NULL)
        display_particles(particles, animation, true,setup->itermax);
}


void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
	}

	// Compute derivatives and normal
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
		// assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
		compute_normal(particles[i], particles_derivatives[i]);
	}

	// Assemble residual and compute curvature
	for (int i = 0; i < n_p; i++) {
	    //particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
	    assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
	}

	// Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
	for (int i = 0; i < n_p; i++)
		time_integrate(particles[i], residuals[i], setup->timestep);
		// time_integrate_CSPM(particles[i],particles_derivatives[i], residuals[i], setup);
}

void compute_Cs(Particle *particle, Kernel kernel, double kh) {
	particle->Cs = 0;
	Particle *pi = particle;
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		particle->Cs += (pj->m / pj->rho) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		node = node->next;
	}
	// printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives) {
	particle->normal = xy_new(0.0, 0.0);
	xy *n = particle_derivatives->grad_Cs; // surface normal inward
	double norm_n = norm(n); // norm of n
	particle->normal->x = n->x / norm_n;
	particle->normal->y = n->y / norm_n;
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual,Setup* setup) {
	double mu_i = particle->param->dynamic_viscosity;

	double rho_i = particle->rho;
	double div_vel_i = particle_derivatives->div_v;
	xy* grad_P = particle_derivatives->grad_P;
	xy* lapl_v = particle_derivatives->lapl_v;

	// Compute UNIT normal vector
	xy *n = particle_derivatives->grad_Cs; // surface normal inward
	double norm_n = norm(n);
	n->x /= norm_n, n->y /= norm_n;

	double lapl_Cs = particle_derivatives->lapl_Cs;
	// Choose between curvature estimated with Laplacian of colour field or with divergence of normal
	// 	double kappa = - lapl_Cs / norm_n; // curvature with Laplacian of colour field

	double fs_x = 0; double fs_y = 0;
	// Apply surface tension only on particles in the vicinity the interface
	bool criterion;
	// Identification based on the norm of the normal
	if (setup->free_surface_detection == CSF)
		criterion = norm_n > setup->interface_threshold;
	// Identification based on the divergence of the position vector
	else if (setup->free_surface_detection == DIVERGENCE)
		criterion = compute_div(particle, Particle_get_pos, setup->kernel, setup->kh) <= setup->interface_threshold;
	else
		criterion = false;
	if (criterion) {
		particle->on_free_surface = true;
	  double kappa = compute_curvature(particle, setup, 0.5);
		particle->P= 0;
	}
	else
		particle->on_free_surface = false;

	residual->mass_eq = -rho_i * div_vel_i;
	residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x;
	residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y;
	if (setup->gravity == 1){
		double g = 9.81;
		residual->momentum_y_eq -= g;
	}
}



// Time integrate the Navier-Stokes equations based on the residual already assembled
void time_integrate(Particle* particle, Residual* residual, double delta_t) {

	// Update density and velocity with an Euler explicit scheme (TODO: implement more accurate and more stable schemes)
	particle->rho += delta_t * residual->mass_eq;
	particle->v->x += delta_t * residual->momentum_x_eq;
	particle->v->y += delta_t * residual->momentum_y_eq;

	// Update position with an Euler Implicit scheme
	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;


	// Update pressure with Tait's equation of state
	double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
	// double B = 0.85*1e5;
	particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);

}
void time_integrate_CSPM(Particle* particle, Particle_derivatives *dp, Residual* residual, Setup* setup) {

	// Update density and velocity with an Euler explicit scheme (TODO: implement more accurate and more stable schemes)
	double delta_t = setup->timestep;
	particle->rho += delta_t * residual->mass_eq;
	Corrective_Smoothed_Particle_Method(particle , dp, setup->kh, setup->kernel);
	particle->v->x += delta_t * residual->momentum_x_eq;
	particle->v->y += delta_t * residual->momentum_y_eq;

	// Update position with an Euler Implicit scheme
	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;


	// Update pressure with Tait's equation of state
	double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
	// double B = 0.85*1e5;
	particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);

}

double compute_curvature(Particle *particle, Setup *setup, double epsilon) {
    double num = 0.0;
    // printf("%lf\n", compute_div(particle, Particle_get_normal, setup->kernel, setup->kh));
    Particle *pi = particle;
    xy *denom = xy_new(0,0);//grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);
    //denom->x=0.0;
    //denom->y=0.0;
    ListNode *node = pi->neighborhood->head;
    double Csi=0.0;

    //Constriction of Csi just like the book
    while (node != NULL) {
        Particle *pj = node->v;
        double mrho= pj->m/pj->rho;
        Csi += mrho*eval_kernel(pi->pos,pj->pos,setup->kh,setup->kernel);
        node = node->next;
    }
    //Constriction of kappa just like the book
    node = pi->neighborhood->head;
    while (node != NULL) {
        Particle *pj = node->v;
        if(pj != pi){
	        //mass density ratio
	        double mrho = pj->m/pj->rho;

	        //nomalized position-> (X_i-X_j) / ||X_i-X_j||^2
	        xy *ij = xy_new(pi->pos->x - pj->pos->x, pi->pos->y - pj->pos->y);
	        double norm_ij = squared(norm(ij));
	        xy *unit = xy_new(ij->x/norm_ij, ij->y/norm_ij);
	        //double unit = (pi->pos-pj->pos)/squared(norm(pi->pos - pj->pos));

	        // Gradient of W
	        xy *grad_W = grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);

	        //Vectorial decompsoition of denom
	        denom->x+=grad_W->x*mrho;
	        denom->y+=grad_W->y*mrho;

	        // Numerator assembly
	        num+=mrho*(Csi-1)*((unit->x*grad_W->x)+(unit->y*grad_W->y));

	        free(grad_W);
    	}
        node = node->next;
    }
    //Kappa assembly
    return (num / norm(denom))/10;
}

void compute_XSPH_correction(Particle *pi, Kernel kernel, double kh, double epsilon) {
	xy_reset(pi->XSPH_correction);
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		pi->XSPH_correction->x += (pj->m / pj->rho) * (pi->v->x - pj->v->x) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		pi->XSPH_correction->y += (pj->m / pj->rho) * (pi->v->y - pj->v->y) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		//printf("%lf\n", pj->m);
		node = node->next;
	}
	pi->XSPH_correction->x *= epsilon;
	pi->XSPH_correction->y *= epsilon;
	//printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

/////////////




double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma) {
  // Relations from "Simulation of surface tension in 2D and 3D with smoothed particle hydrodynamics method", Zhang (2010)
  double dt_1 = 0.25 * h_p / c_0; // propagation of sound waves
  double dt_2 = INFINITY;
  if (mu > 0.0) dt_2 = 0.25 * (h_p*h_p) / (mu / rho_0); // viscous diffusion
  double dt_3 = INFINITY;
  if (sigma > 0.0) dt_3 = 0.25 * sqrt((rho_0*pow(h_p,3))/(2*M_PI*sigma)); // surface tension (capillary waves)

  double dt_min_interm = fmin(dt_1, dt_2);
  return safety_param * fmin(dt_min_interm, dt_3);

}

Boundary* Boundary_new(double xleft1, double xleft2,  double xright, double ybottom1, double ybottom2, double ytop,double CR, double CF){
	Boundary* boundary = (Boundary*) malloc(sizeof(Boundary));
	boundary->xleft1 = xleft1;
	boundary->xleft2 = xleft2;
	boundary->xright = xright;
	boundary->ybottom1 = ybottom1;
	boundary->ybottom2 = ybottom2;
	boundary->ytop = ytop;
	boundary->CR = CR;
	boundary->CF = CF;
	return boundary;
}
void Boundary_free(Boundary* boundary){
	free(boundary);
}

void center_reflection_right(Particle* pi, double CR, double Rp, double d){
		pi->pos->x -= (1+CR)*(Rp - d);
}
void center_reflection_left(Particle* pi, double CR, double Rp, double d){
		pi->pos->x += (1+CR)*(Rp - d);
}
void center_reflection_top(Particle* pi, double CR, double Rp, double d){
		pi->pos->y -= (1+CR)*(Rp - d);
}
void center_reflection_bottom(Particle* pi, double CR, double Rp, double d){
		pi->pos->y += (1+CR)*(Rp - d);
}
void velocity_reflection_vertical(Particle* pi, double CR, double CF){
	double vpN = pi->v->x;
	double vpT = pi->v->y;
	pi->v->x = -vpN*CR;
	pi->v->y = (1-CF)* vpT;
}
void velocity_reflection_horizontal(Particle* pi, double CR, double CF){
	double vpN = pi->v->y;
	double vpT = pi->v->x;
	pi->v->y = -vpN*CR;
	pi->v->x = (1-CF)* vpT;
}
void reflective_boundary(Particle** p, int n_p, Boundary* boundary, double Rp){
	// We have just computed the time integration. We correct the positions of the particles
	printf("je suis rentre dans reflective boundary \n");
	double CF = boundary->CF;
	double CR = boundary->CR;
	for(int i = 0; i < n_p ; i++){
		// For each particle we check if its position is close to a wall
		Particle* pi = p[i];
		double xright = boundary->xright;	double xleft1 = boundary->xleft1;	double xleft2 = boundary->xleft2;
		double ytop = boundary->ytop;			double ybottom1 = boundary->ybottom1;		double ybottom2 = boundary->ybottom2;
		double dxright,dxleft,dytop,dybottom;

		// Collision test
		int collision, collision_right, collision_left, collision_top, collision_bottom;
		collision = 0; collision_right = 0; collision_left = 0; collision_top = 0; collision_bottom = 0;

		int i = 0;
		while(pi->pos->x > boundary->xright - Rp){
			collision_right = 1;
			collision = 1;
			dxright = fabs(pi->pos->x - xright);
			center_reflection_right(pi,CR,Rp,dxright);
			velocity_reflection_vertical(pi,CR,CF);
			i++;
			// printf("i = %d\n",i);
			// printf("position : %f,%f\n", pi->pos->x, pi->pos->y);
		}
		while(pi->pos->x < boundary->xleft1 + Rp && pi->pos->y > boundary->ybottom1 + Rp){
			collision_left = 1;
			collision = 1;
			dxleft = fabs(pi->pos->x - xleft1);
			center_reflection_left(pi,CR,Rp,dxleft);
			velocity_reflection_vertical(pi,CR,CF);

		}
		while (pi->pos->x < boundary->xleft2 + Rp && pi->pos->y <= boundary->ybottom1 ) {
			collision_left = 1;
			collision = 1;
			dxleft = fabs(pi->pos->x - xleft2);
			center_reflection_left(pi, CR, Rp, dxleft);
			velocity_reflection_vertical(pi, CR, CF);

		}
		while(pi->pos->y > boundary->ytop - Rp){
			collision_top = 1;
			collision = 1;
			dytop = fabs(pi->pos->y - ytop);
			center_reflection_top(pi,CR,Rp,dytop);
			velocity_reflection_horizontal(pi,CR,CF);
		}
		while(pi->pos->y < boundary->ybottom1 + Rp && pi->pos->x < boundary->xleft2){
			collision_bottom = 1;
			collision = 1;
			dybottom = fabs(pi->pos->y - ybottom1);
			center_reflection_bottom(pi,CR,Rp,dybottom);
			velocity_reflection_horizontal(pi,CR,CF);
		}
		while (pi->pos->y < boundary->ybottom2 + Rp && pi->pos->x >= boundary->xleft2) {
			collision_bottom = 1;
			collision = 1;
			dybottom = fabs(pi->pos->y - ybottom2);
			center_reflection_bottom(pi, CR, Rp, dybottom);
			velocity_reflection_horizontal(pi, CR, CF);
		}
	}
}
