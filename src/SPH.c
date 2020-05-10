#include "SPH.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* These are the main functions of the code. */




/* Main for loop on the time iteration. The steps for the computation are the following :
 1) update of the cells : which particles are in which cells
 2) update of the neighborhoods of each particle
 3) update of the position of each particle based on the resolution of the Navier - Stokes equations
 4) Check for boundary treatment through the reflective boundary approach.
 */
void simulate_SPH(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup, Animation* animation, Boundary* boundary){
    double current_time = 0.0;
    double Rp = 0.001;
    int ii = 5;
    double bounds[6] = { boundary->xleft1,boundary->xleft2, boundary->xright, boundary->ybottom1, boundary->ybottom2, boundary->ytop };
    printf("%d\n", setup->itermax);
    
    
    FILE* myfile = fopen("Energy.txt", "w");
    if (myfile == NULL)
    {
        printf("Error opening file!\n");
    }
    double t_neighborhood = 0;
    for (int iter = 0; iter < setup->itermax; iter++) {
        printf("----------------------------------------------------- \n");
        printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
        if (iter%setup->verlet->T == 0) {
            double v_max = max_velocity(particles, n_p);
            printf("v_max = %f \n",v_max);
            setup->verlet->L = 1.1*(2*v_max*(double)setup->verlet->T*setup->timestep);;
        }
        clock_t t0 = clock();
        if (setup->verlet->use_verlet) {
            update_verlet_cells(grid, particles, n_p, setup->verlet);
        }
        else {
            update_cells(grid, particles, n_p);
        }
        update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
        clock_t t1 = clock();
        t_neighborhood += (double)(t1-t0);
        if (animation != NULL ) { //&& iter%1 == 0)
            double Epot = 0.;
            double Ekin = 0.;
            double vel;
            for (int i=0; i<n_p; i++) {
                Epot += (particles[i]->pos->y+0.5)*particles[i]->m*9.81;
                vel = particles[i]->v->x*particles[i]->v->x + particles[i]->v->y*particles[i]->v->y;
                Ekin += vel*particles[i]->m / 2.;
            }
            fprintf(myfile,"%0.12f %0.12f\n",Epot, Ekin);
            fflush(myfile);
            display_particles_boundary(particles, animation, false,iter,bounds);
             
        }
        update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);
        reflective_boundary(particles,n_p,boundary,Rp);
        current_time += setup->timestep;
    }
    printf("total neighborhood search time : %f\n",t_neighborhood/CLOCKS_PER_SEC);
    fclose(myfile);
}

/* Position update of the particles in different steps :
 1) Computation of the color field of the particles
 2) computation of the spatial derivatives
 3) computation of the momentum and mass balance equations
 4) time integration and application of the xpsh correction
 */
void update_positions(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {
    
    for (int i = 0; i < n_p; i++) {
        compute_Cs(particles[i], setup->kernel, setup->kh);
    }
    
    for (int i = 0; i < n_p; i++) {
        compute_derivatives(particles[i], particles_derivatives[i], setup->kernel, setup->kh);
        compute_normal(particles[i], particles_derivatives[i]);
    }
    
    for (int i = 0; i < n_p; i++) {
        assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
    }
    
    for (int i = 0; i < n_p; i++) {
        compute_XSPH_correction(particles[i], setup->kernel, setup->kh, 0.5);
        time_integrate(particles[i], residuals[i], setup->timestep);
    }
}

/*
 Application of the balance equations in momentum and mass. Also verification if the particles are on a free surface or not.
 */
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual,Setup* setup) {
    double mu_i = particle->param->dynamic_viscosity;

    double rho_i = particle->rho;
    double div_vel_i = particle_derivatives->div_v;
    xy* grad_P = particle_derivatives->grad_P;
    xy* lapl_v = particle_derivatives->lapl_v;
    xy* art_visc = particle_derivatives->art_visc;

    xy *n = particle_derivatives->grad_Cs;
    double norm_n = norm(n);
    n->x /= norm_n, n->y /= norm_n;

    double lapl_Cs = particle_derivatives->lapl_Cs;

    double fs_x = 0; double fs_y = 0;
    
    bool criterion = compute_div(particle, Particle_get_pos, setup->kernel, setup->kh) <= setup->interface_threshold;
    
    if (criterion && !particle->on_boundary) {
        particle->on_free_surface = true;
        xy* fs = compute_surfaceTension(particle, particle_derivatives);
        particle->P = 0;
        residual->mass_eq = -rho_i * div_vel_i;
        residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x - art_visc->x + fs->x;
        residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y - art_visc->y + fs->y;
        free(fs);
    }
    else {
        particle->on_free_surface = false;
        residual->mass_eq = -rho_i * div_vel_i;
        residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x - art_visc->x;
        residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y - art_visc->y;
    }
    if (setup->gravity == 1){
        double g = 9.81;
        residual->momentum_y_eq -= g;
    }
}

/*
 Time integration of the Navier-Stokes equations based on the residuals assembled earlier.
 The scheme is an Euler explicit scheme.
 */
void time_integrate(Particle* particle, Residual* residual, double delta_t) {

    particle->rho += delta_t * residual->mass_eq;
    particle->v->x += delta_t * residual->momentum_x_eq;
    particle->v->y += delta_t * residual->momentum_y_eq;

    particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
    particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;
    
    double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
    if (particle->on_free_surface && !particle->on_boundary) {
        particle->P = 0.;
    }
    else {
        particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);
    }
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
    Node_free(node);
}

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives) {
    xy *n = particle_derivatives->grad_Cs;
    particle->normal->x = n->x;
    particle->normal->y = n->y;
}

xy* compute_surfaceTension(Particle* p, Particle_derivatives* d)
{
    double sigma = 72.86e-3;
    double lapl_Cs = d->lapl_Cs;
    double n_x = p->normal->x;
    double n_y = p->normal->y;
    double norm = sqrt(n_x*n_x + n_y*n_y);
    double fs_x = -sigma*lapl_Cs*n_x/norm;
    double fs_y = -sigma*lapl_Cs*n_y/norm;
    if (norm>20.) {
        return xy_new(fs_x,fs_y);
    }
    else{
        p->on_free_surface = false;
        return xy_new(0.,0.);
    }
}

void compute_XSPH_correction(Particle *pi, Kernel kernel, double kh, double epsilon) {
    xy_reset(pi->XSPH_correction);
    ListNode *node = pi->neighborhood->head;
    while (node != NULL) {
        Particle *pj = node->v;
        pi->XSPH_correction->x += (2*pj->m) / (pj->rho+pi->rho) * (pj->v->x - pi->v->x) * eval_kernel(pi->pos, pj->pos, kh, kernel);
        pi->XSPH_correction->y += (2*pj->m) / (pj->rho+pi->rho) * (pj->v->y - pi->v->y) * eval_kernel(pi->pos, pj->pos, kh, kernel);
        node = node->next;
    }
    Node_free(node);
    pi->XSPH_correction->x *= epsilon;
    pi->XSPH_correction->y *= epsilon;
}


void center_reflection_right(Particle* pi, double CR, double Rp, double d){
        pi->pos->x -= (1+CR)*(d);
}
void center_reflection_left(Particle* pi, double CR, double Rp, double d){
        pi->pos->x += (1+CR)*(d);
}
void center_reflection_top(Particle* pi, double CR, double Rp, double d){
        pi->pos->y -= (1+CR)*(d);
}
void center_reflection_bottom(Particle* pi, double CR, double Rp, double d){
        pi->pos->y += (1+CR)*(d);
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
    double CF = boundary->CF;
    double CR = boundary->CR;
    double xright = boundary->xright;    double xleft1 = boundary->xleft1;    double xleft2 = boundary->xleft2;
    double ytop = boundary->ytop;            double ybottom1 = boundary->ybottom1;        double ybottom2 = boundary->ybottom2;
    double dxright,dxleft,dytop,dybottom;
    double tol = fabs(0.01*(xleft1-xright));
    for(int i = 0; i < n_p ; i++){
        Particle* pi = p[i];
        pi->on_boundary=false;
        double x = pi->pos->x;
        double y = pi->pos->y;
        if ( fabs(x-xright) < tol || (fabs(y-ybottom2) < tol && x > xleft2) || (fabs(y-ybottom1) < tol && x < xleft2)) {
            pi->on_boundary = true;
        }
        bool bounced = false;
        while (!bounced) {
            bounced = true;
            if (pi->pos->x > boundary->xright && bounced) {
                dxright = fabs(pi->pos->x - xright);
                center_reflection_right(pi, CR, Rp, dxright);
                velocity_reflection_vertical(pi, CR, CF);
                pi->on_boundary=true;
                bounced = false;
            }
            if (pi->pos->x < boundary->xleft1 && bounced) {
                dxleft = fabs(pi->pos->x - xleft1);
                center_reflection_left(pi, CR, Rp, dxleft);
                velocity_reflection_vertical(pi, CR, CF);
                pi->on_boundary=true;
                bounced = false;
            }
            if (pi->pos->y > boundary->ytop && bounced) {
                dytop = fabs(pi->pos->y - ytop);
                center_reflection_top(pi, CR, Rp, dytop);
                velocity_reflection_horizontal(pi, CR, CF);
                pi->on_boundary=true;
                bounced = false;
            }
            if (pi->pos->y < boundary->ybottom2 && bounced) {
                dybottom = fabs(pi->pos->y - ybottom2);
                center_reflection_bottom(pi, CR, Rp, dybottom);
                velocity_reflection_horizontal(pi, CR, CF);
                pi->on_boundary=true;
                bounced = false;
            }
            if (pi->pos->x < boundary->xleft2 && pi->pos->y < boundary->ybottom1 && bounced) {
                if (fabs(pi->pos->x-boundary->xleft2) >= fabs(pi->pos->y - boundary->ybottom1)) {
                    dybottom = fabs(pi->pos->y - ybottom1);
                    center_reflection_bottom(pi, CR, Rp, dybottom);
                    velocity_reflection_horizontal(pi, CR, CF);
                    pi->on_boundary=true;
                }
                else if (fabs(pi->pos->y - boundary->ybottom1) > fabs(pi->pos->x - boundary->xleft2)) {
                    dxleft = fabs(pi->pos->x - xleft2);

                    center_reflection_left(pi, CR, Rp, dxleft);
                    velocity_reflection_vertical(pi, CR, CF);
                    pi->on_boundary=true;
                }
                else {
                    dybottom = fabs(pi->pos->y - ybottom1);
                    center_reflection_bottom(pi, CR, Rp, dybottom);
                    velocity_reflection_horizontal(pi, CR, CF);
                    pi->on_boundary=true;
                }
                bounced = false;
            }
        }
    }
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
    double pos_x;
    double pos_y;
    int index;
    for(int i = 0 ; i < n_p; i++){
        double v = sqrt(p[i]->v->x*p[i]->v->x + p[i]->v->y*p[i]->v->y);
        if(max < v){
            max = v;
            pos_x = p[i]->pos->x;
            pos_y = p[i]->pos->y;
            index = p[i]->index;
        }
    }
    printf("max velocity particle : %d, position : (%f,%f) \n",index,pos_x,pos_y);
    return max;
}
