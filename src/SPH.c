#include "SPH.h"

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
        double v = sqrt(p[i]->v->x*p[i]->v->x + p[i]->v->y*p[i]->v->y);
        if(max < v){
            max = v;
        }
    }
    return max;
}

void simulate_boundary(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation, Boundary* boundary){
    double current_time = 0.0;
    double Rp = 0.001;
    int ii = 5;
    double bounds[6] = { boundary->xleft1,boundary->xleft2, boundary->xright, boundary->ybottom1, boundary->ybottom2, boundary->ytop };
    printf("%d\n", setup->itermax);
    for (int iter = 0; iter < setup->itermax; iter++) {
        printf("----------------------------------------------------- \n");
        printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
        if (iter%setup->verlet->T == 0) {
            double v_max = max_velocity(particles, n_p);
            printf("v_max = %f \n",v_max);
            setup->verlet->L = 1.3*v_max;
        }
        update_cells(grid, particles, n_p);
        update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
        if (animation != NULL && iter%1 == 0)
            display_particles_boundary(particles, animation, false,iter,bounds);
        update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);

        //printf("velocity_max = %f\n", max_velocity(particles,n_p));
        reflective_boundary(particles,n_p,boundary,Rp);
        if (iter%ii == 0){
            // density_correction_MLS(particles, n_p, setup->kh, setup->kernel);
        }
        current_time += setup->timestep;
    }
}

void simulate_boundary_flow(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation, Boundary* boundary) {
    double current_time = 0.0;
    double Rp = 0.001;
    int ii = 5;
    double bounds[6] = { boundary->xleft1,boundary->xleft2, boundary->xright, boundary->ybottom1, boundary->ybottom2, boundary->ytop };
    printf("%d\n", setup->itermax);
    for (int iter = 0; iter < setup->itermax; iter++) {
        printf("----------------------------------------------------- \n");
        printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
        if (iter % 40 == 0 && iter < 40*20) {
            double l = 0.057;
            double h = 2 * l / 30;
            for (int num = 0; num < 30; num++) {
                int index1 = n_p + num;
                xy* pos1;
                xy* v1;
                pos1 = xy_new(-2 * l + 0.001, num * h);
                v1 = xy_new(2, 0.0); // initial velocity = 0
                particles[index1]->pos = pos1;
                particles[index1]->v = v1;
            }
            n_p += 30;
        }
        update_cells(grid, particles, n_p);
        update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
        //printf("j'ai update neighboor \n");
        if (animation != NULL)
            display_particles_boundary(particles, animation, false, iter, bounds);
        update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);

        //printf("velocity_max = %f\n", max_velocity(particles, n_p));
        reflective_boundary(particles, n_p, boundary, Rp);
        if (iter % ii == 0) {
            // density_correction_MLS(particles, n_p, setup->kh, setup->kernel);
        }
        //get_M0(particles,n_p,setup->kh,setup->kernel);
        //get_M1(particles,n_p,setup->kh,setup->kernel);
        current_time += setup->timestep;
    }
}



void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

    // Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
    //for (int i = 0; i < n_p; i++) {
        //compute_Cs(particles[i], setup->kernel, setup->kh);
        //if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
    //}
    for (int i = 0; i < n_p; i++) {
        compute_Cs(particles[i], setup->kernel, setup->kh);
    }
    // Compute derivatives and normal
    for (int i = 0; i < n_p; i++) {
        compute_derivatives(particles[i], particles_derivatives[i], setup->kernel, setup->kh);
        compute_normal(particles[i], particles_derivatives[i]);
    }
    
    for (int i = 0; i < n_p; i++) {
        assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
    }
    
    // Apply XSPH correction and time integrate
    for (int i = 0; i < n_p; i++) {
        compute_XSPH_correction(particles[i], setup->kernel, setup->kh, 0.5);
        time_integrate(particles[i], residuals[i], setup->timestep);
    }
     
    // Assemble residual and compute curvature
    /*
     
    
    for (int i = 0; i < n_p; i++) {
        //particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
        
    }

    // Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
    for (int i = 0; i < n_p; i++)
        // time_integrate_CSPM(particles[i],particles_derivatives[i], residuals[i], setup);
     */
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
    // printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives) {
    //particle->normal = xy_new(0.0, 0.0);
    xy *n = particle_derivatives->grad_Cs; // surface normal inward
    //double norm_n = norm(n); // norm of n
    particle->normal->x = n->x;
    particle->normal->y = n->y;
}

xy* compute_surfaceTension(Particle* p, Particle_derivatives* d)
{
    double sigma = 72.86e-3;
    double lapl_Cs = d->lapl_Cs;
    //printf("laplacianCs = %f\n",lapl_Cs);
    double n_x = p->normal->x;
    double n_y = p->normal->y;
    double norm = sqrt(n_x*n_x + n_y*n_y);
    //printf("normal norm : %f\n",norm);
    double fs_x = -sigma*lapl_Cs*n_x/norm;
    double fs_y = -sigma*lapl_Cs*n_y/norm;
    //printf("(%f,%f)\n",fs_x,fs_y);
    if (norm>20.) {
        return xy_new(fs_x,fs_y);
    }
    else{
        return xy_new(0.,0.);
    }
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual,Setup* setup) {
    double mu_i = particle->param->dynamic_viscosity;

    double rho_i = particle->rho;
    double div_vel_i = particle_derivatives->div_v;
    xy* grad_P = particle_derivatives->grad_P;
    xy* lapl_v = particle_derivatives->lapl_v;
    xy* art_visc = particle_derivatives->art_visc;

    // Compute UNIT normal vector
    xy *n = particle_derivatives->grad_Cs; // surface normal inward
    double norm_n = norm(n);
    n->x /= norm_n, n->y /= norm_n;

    double lapl_Cs = particle_derivatives->lapl_Cs;
    // Choose between curvature estimated with Laplacian of colour field or with divergence of normal
    //     double kappa = - lapl_Cs / norm_n; // curvature with Laplacian of colour field

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
    if (criterion && particle->on_boundary==false) {
        particle->on_free_surface = true;
        //double kappa = compute_curvature(particle, setup, 0.5);
        xy* fs = compute_surfaceTension(particle, particle_derivatives);
        particle->P = 0;
        residual->mass_eq = -rho_i * div_vel_i;
        residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x + fs->x - art_visc->x;
        residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y + fs->y - art_visc->y;
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
    //double B = 0.85*1e5;
    if (particle->on_free_surface && !particle->on_boundary) {
        particle->P = 0.;
    }
    else {
        particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1);
    }
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
    Particle *pi = particle;
    xy *denom = xy_new(0,0);//grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);
    //denom->x=0.0;
    //denom->y=0.0;
    ListNode *node = pi->neighborhood->head;
    double Csi=0.0;

    //Constriction of Csi just like the book
    /*
    while (node != NULL) {
        Particle *pj = node->v;
        double mrho= pj->m/pj->rho;
        Csi += mrho*eval_kernel(pi->pos,pj->pos,setup->kh,setup->kernel);
        node = node->next;
    }
     */
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
            free(ij);
            free(unit);
        }
        node = node->next;
    }
    Node_free(node);
    free(denom);
    //Kappa assembly
    return (num / norm(denom))/10;
}

void compute_XSPH_correction(Particle *pi, Kernel kernel, double kh, double epsilon) {
    xy_reset(pi->XSPH_correction);
    ListNode *node = pi->neighborhood->head;
    while (node != NULL) {
        Particle *pj = node->v;
        pi->XSPH_correction->x += (2*pj->m) / (pj->rho+pi->rho) * (pj->v->x - pi->v->x) * eval_kernel(pi->pos, pj->pos, kh, kernel);
        pi->XSPH_correction->y += (2*pj->m) / (pj->rho+pi->rho) * (pj->v->y - pi->v->y) * eval_kernel(pi->pos, pj->pos, kh, kernel);
        //printf("%lf\n", pj->m);
        node = node->next;
    }
    Node_free(node);
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
    //printf("reflective boundary begin \n");
    // We have just computed the time integration. We correct the positions of the particles
    double CF = boundary->CF;
    double CR = boundary->CR;
    double xright = boundary->xright;    double xleft1 = boundary->xleft1;    double xleft2 = boundary->xleft2;
    double ytop = boundary->ytop;            double ybottom1 = boundary->ybottom1;        double ybottom2 = boundary->ybottom2;
    double dxright,dxleft,dytop,dybottom;
    for(int i = 0; i < n_p ; i++){
        // For each particle we check if its position is close to a wall
        Particle* pi = p[i];
        pi->on_boundary=false;

        bool bounced = false;
        //while (!bounced) {
            bounced = true;
            if (pi->pos->x > boundary->xright) {
                //printf("in 1 \n");
                dxright = fabs(pi->pos->x - xright);
                center_reflection_right(pi, CR, Rp, dxright);
                velocity_reflection_vertical(pi, CR, CF);
                pi->on_boundary=true;
                //printf("out 1 \n");
                bounced = false;
            }
            if (pi->pos->x < boundary->xleft1 ) {
                //printf("in 2 \n");
                dxleft = fabs(pi->pos->x - xleft1);
                center_reflection_left(pi, CR, Rp, dxleft);
                velocity_reflection_vertical(pi, CR, CF);
                pi->on_boundary=true;
                //printf("out 2 \n");
                bounced = false;
            }
            if (pi->pos->y > boundary->ytop) {
                //printf("in 4 \n");
                dytop = fabs(pi->pos->y - ytop);
                center_reflection_top(pi, CR, Rp, dytop);
                velocity_reflection_horizontal(pi, CR, CF);
                pi->on_boundary=true;
                //printf("out 4 \n");
                bounced = false;
            }
            if (pi->pos->y < boundary->ybottom2) {
                //printf("in 6 \n");
                dybottom = fabs(pi->pos->y - ybottom2);
                center_reflection_bottom(pi, CR, Rp, dybottom);
                velocity_reflection_horizontal(pi, CR, CF);
                pi->on_boundary=true;
                //printf("out 6 \n");
                bounced = false;
            }
            if (pi->pos->x < boundary->xleft2 && pi->pos->y < boundary->ybottom1 ) {
                if (fabs(pi->pos->x-boundary->xleft2) >= fabs(pi->pos->y - boundary->ybottom1)) {
                    //printf("in 5 \n");
                    //printf("x : %f\n", pi->pos->x);
                    //printf("y : %f\n", pi->pos->y);
                    dybottom = fabs(pi->pos->y - ybottom1);
                    center_reflection_bottom(pi, CR, Rp, dybottom);
                    velocity_reflection_horizontal(pi, CR, CF);
                    pi->on_boundary=true;
                    //printf("out 5 \n");
                }
                else if (fabs(pi->pos->y - boundary->ybottom1) > fabs(pi->pos->x - boundary->xleft2)) {
                    //printf("in 3 \n");
                    //printf("ind pi : %d\n", pi->index);
                    //printf("initial x : %f\n", pi->pos->x);
                    //printf("xleft2 = %f \n", xleft2);
                    dxleft = fabs(pi->pos->x - xleft2);

                    center_reflection_left(pi, CR, Rp, dxleft);
                    velocity_reflection_vertical(pi, CR, CF);
                    pi->on_boundary=true;
                    //printf("dxleft : %f\n", dxleft);
                    //printf("x : %f\n", pi->pos->x);
                    //printf("y : %f\n", pi->pos->y);
                    //printf("out 3 \n");
                }
                else {
                    //printf("in 5 but too much in the middle -> reduce dt\n");
                    dybottom = fabs(pi->pos->y - ybottom1);
                    center_reflection_bottom(pi, CR, Rp, dybottom);
                    velocity_reflection_horizontal(pi, CR, CF);
                    pi->on_boundary=true;
                    //printf("out 5 \n");
                }
                bounced = false;
            }
    }
}
