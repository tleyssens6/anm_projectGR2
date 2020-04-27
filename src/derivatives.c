#include "derivatives.h"

void compute_derivatives(Particle* p, Particle_derivatives* deriv, Kernel kernel, double kh)
{
    ListNode *node = p->neighborhood->head;
    deriv->div_v = 0.;
    xy_reset(deriv->lapl_v);
    xy_reset(deriv->grad_P);
    xy_reset(deriv->grad_Cs);
    deriv->lapl_Cs = 0.;
    xy *grad_W = xy_new(0.,0.);
    xy *DXij = xy_new(0.,0.);
    xy* v_q;
    xy* v_p = Particle_get_v(p);
    double P_p = Particle_get_P(p);
    double P_q;
    
    double Cs_p = Particle_get_Cs(p);
    double Cs_q;
    
    double d2;
    
    while(node != NULL) {
        Particle *q = node->v;
        xy *grad_W = grad_kernel(p->pos, q->pos, kh, kernel);
        v_q = Particle_get_v(q);
        deriv->div_v += ((v_p->x - v_q->x) * grad_W->x + (v_p->y - v_q->y) * grad_W->y) * q->m/p->rho;
        d2 = squared(p->pos->x - q->pos->x) + squared(p->pos->y - q->pos->y);
        xy *DXij = xy_new((p->pos->x - q->pos->x) / d2, (p->pos->y - q->pos->y) / d2);
        deriv->lapl_v->x += 2 * (q->m/q->rho) * (v_p->x - v_q->x) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        deriv->lapl_v->y += 2 * (q->m/q->rho) * (v_p->y - v_q->y) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        P_q = Particle_get_P(q);
        deriv->grad_P->x -= p->rho * q->m * (P_p/squared(p->rho) + P_q/squared(q->rho)) * grad_W->x;
        deriv->grad_P->y -= p->rho * q->m * (P_p/squared(p->rho) + P_q/squared(q->rho)) * grad_W->y;
        Cs_q = Particle_get_Cs(q);
        deriv->grad_Cs->x -= p->rho * q->m * (Cs_p/squared(p->rho) + Cs_q/squared(q->rho)) * grad_W->x;
        deriv->grad_Cs->y -= p->rho * q->m * (Cs_p/squared(p->rho) + Cs_q/squared(q->rho)) * grad_W->y;
        deriv->lapl_Cs += 2 * (q->m/q->rho) * (Cs_p - Cs_q) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        node = node->next;
        free(grad_W);
        free(DXij);
    }
    //free(grad_W);
    //free(DXij);
    //free(v_q);
    //free(v_p);
}



double compute_div(Particle* pi, xy_getter get, Kernel kernel, double kh) {
    double div = 0;
    xy *fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;

        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        // correct_grad_local(grad_W, pi, pj, kh, kernel);

        xy *fj = get(pj);
        div += ((fj->x - fi->x) * grad_W->x)*pj->m;
        div += ((fj->y - fi->y) * grad_W->y)*pj->m;
        free(grad_W);
        node = node->next;
    }
    div /= pi->rho;
    return div;
}

void compute_grad(Particle* pi, scalar_getter get, Kernel kernel, double kh,xy* grad) {

	double gx = 0;
	double gy = 0;
  double fi = get(pi);
  ListNode *node = pi->neighborhood->head;
  //printf("Computing gradient of (%lf, %lf), fi = %lf\n", particle->pos->x, particle->pos->y, fi);
  while(node != NULL) {
      Particle *pj = node->v;
      xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);

      double fj = get(pj);
      //printf("Position of pj: (%lf, %lf), fj = %lf\n", pj->pos->x, pj->pos->y, fj);
      gx += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->x; // sign is not the same as in the def...
      gy += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->y;
      free(grad_W);
      //printf("grad = (%lf, %lf), fj = %lf\n", grad->x, grad->y, fj);
      node = node->next;
  }
  // xy *cg = xy_new(gx,gy);
  grad->x = gx;
  grad->y = gy;
  // correct_grad(grad, pi, kh, kernel);

  // free(cg);
}

double compute_lapl(Particle *pi, scalar_getter get, Kernel kernel, double kh) {
    double lapl = 0;
    double fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        // correct_grad_local(grad_W, pi, pj, kh, kernel);

        double fj = get(pj);
        double d2 = squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y); // squared distance between particles
        // xy *DXij = xy_new((pi->pos->x - pj->pos->x) / d2, (pi->pos->y - pj->pos->y) / d2); // Delta X_{ij}
        xy *DXij = xy_new((pi->pos->x - pj->pos->x) / d2, (pi->pos->y - pj->pos->y) / d2); // WARNING
        if(d2 != 0) lapl += 2 * (pj->m/pj->rho) * (fi - fj) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        free(grad_W);
        free(DXij);
        node = node->next;
    }
    return lapl;
}
