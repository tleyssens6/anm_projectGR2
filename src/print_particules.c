#include "print_particules.h"

void colormap_Cs(Particle *p, float color[3]);

// Fills data with particle data
void fillData(GLfloat(*data_pool)[8], GLfloat(*data_flow)[8], Particle** particles, int N_pool, int N_flow) {
// 	float rmax = 100.0*sqrtf(2.0f);
    /*
	double max_norm_fs = 0.0;
	double fs_norm_local;
	for (int i = 0; i < N_pool; i++) {
	    xy* fs = xy_new(-particles[i]->param->sigma * particles[i]->normal->x * particles[i]->kappa, -particles[i]->param->sigma * particles[i]->normal->y * particles[i]->kappa);
	    fs_norm_local = norm(fs);
	    if (fs_norm_local > max_norm_fs) max_norm_fs = fs_norm_local;
	}
    for (int i = N_pool; i < N_flow; i++) {
        xy* fs = xy_new(-particles[i]->param->sigma * particles[i]->normal->x * particles[i]->kappa, -particles[i]->param->sigma * particles[i]->normal->y * particles[i]->kappa);
        fs_norm_local = norm(fs);
        if (fs_norm_local > max_norm_fs) max_norm_fs = fs_norm_local;
    }
    */
	for (int i = 0; i < N_pool; i++) {
		Particle* p = particles[i];
		data_pool[i][0] = p->pos->x;
		data_pool[i][1] = p->pos->y;
		data_pool[i][2] = p->v->x;
		data_pool[i][3] = p->v->y;
		// colormap_cell(p, &data[i][4]); // fill color
		//colormap_Cs(p, &data[i][4]); // fill color
        if (p->on_free_surface) {
            //           colormap_uni_color_2(&data[i][4]);
            //colormap_fs(p, &data[i][4], max_norm_fs);
            data_pool[i][4] = 1.0;
            data_pool[i][5] = 1.0;
            data_pool[i][6] = 0.0;
        }
        else {
            //colormap_uni_color(&data[i][4]);
            data_pool[i][4] = 0.6;
            data_pool[i][5] = 0.0;
            data_pool[i][6] = 0.0;
        }

		data_pool[i][7] = 0.8f; // transparency
	}
    for (int i = 0; i < N_flow; i++) {
        Particle* p = particles[i+N_pool];
        data_flow[i][0] = p->pos->x;
        data_flow[i][1] = p->pos->y;
        data_flow[i][2] = p->v->x;
        data_flow[i][3] = p->v->y;
        // colormap_cell(p, &data[i][4]); // fill color
        //colormap_Cs(p, &data[i][4]); // fill color
        if (p->on_free_surface) {
//           colormap_uni_color_2(&data[i][4]);
          //colormap_fs(p, &data[i][4], max_norm_fs);
            data_flow[i][4] = 0.0;
            data_flow[i][5] = 1.0;
            data_flow[i][6] = 1.0;
        }
        else {
          //colormap_uni_color(&data_flow[i][4]);
            data_flow[i][4] = 0.0;
            data_flow[i][5] = 0.0;
            data_flow[i][6] = 0.6;
        }

        data_flow[i][7] = 0.8f; // transparency
    }
}

bov_points_t * load_Grid(Grid* grid,double scale)
{
	int nLines = (grid->nCellx + 1) + (grid->nCelly + 1);
	printf("\n%d\n", nLines);
	GLfloat(*data)[2] = malloc(sizeof(data[0]) * 2 * nLines);
	for (int i = 0;i < (grid->nCellx + 1);i++)
	{
		data[2 * i][0] = grid->left + i * grid->h;
		data[2 * i][1] = grid->bottom;
		data[2 * i + 1][0] = grid->left + i * grid->h;
		data[2 * i + 1][1] = grid->top;
	}
	int shift = 2 * (grid->nCellx + 1);
	for (int j = 0;j < (grid->nCelly + 1);j++)
	{
		data[shift + 2 * j][0] = grid->left;
		data[shift + 2 * j][1] = grid->bottom + j * grid->h;
		data[shift + 2 * j + 1][0] = grid->right;
		data[shift + 2 * j + 1][1] = grid->bottom + j * grid->h;
	}
	bov_points_t *points = bov_points_new(data, 2 * nLines, GL_STATIC_DRAW);
	bov_points_set_width(points, 0.0005);
	double L = grid->h*grid->nCellx;
	bov_points_scale(points, (GLfloat[2]){0.8/L*scale, 0.8/L*scale});
	//bov_points_scale(points, (GLfloat[2]) { 0.008, 0.008 });
	free(data);
	return points;
}

Animation* Animation_new(int N_pool, int N_flow, double timeout,Grid* grid,double scale)
{
	Animation* animation = (Animation*)malloc(sizeof(Animation));
	animation->window = bov_window_new(1024, 780, "ANM Project: SPH");
    bov_window_set_color(animation->window, (GLfloat[]) { 1.f, 1.f, 1.f, 0.0f });
	bov_window_enable_help(animation->window);
	animation->N_pool = N_pool;
    animation->N_flow = N_flow;
	animation->timeout = timeout;
	double L = grid->h*grid->nCellx;
	////set-up particles////
    
	GLfloat(*data_flow)[8] = malloc(sizeof(data_flow[0])*N_flow);
    GLfloat(*data_pool)[8] = malloc(sizeof(data_pool[0])*N_pool);
	bov_points_t *particles_pool = bov_particles_new(data_pool, N_pool, GL_STATIC_DRAW);
    bov_points_t *particles_flow = bov_particles_new(data_flow, N_flow, GL_STATIC_DRAW);
	free(data_flow);
    free(data_pool);
    
	// setting particles appearance
	bov_points_set_width(particles_pool, 0.005);
    bov_points_set_width(particles_flow, 0.005);
	bov_points_set_outline_width(particles_pool, 0.00025);
    bov_points_set_outline_width(particles_flow, 0.00025);

	double c = 4;
	// bov_points_scale(particles, (GLfloat[2]){0.8*c/L*scale, 0.8*c/L*scale});//0.8
	//bov_points_scale(particles, (GLfloat[2]){ 0.008, 0.008 });
	animation->particles_pool = particles_pool;
    animation->particles_flow = particles_flow;
	////set-up grid////
	if (grid != NULL)
		animation->grid = load_Grid(grid,scale*c);
	else
		animation->grid = NULL;
	return animation;
}
void Animation_free(Animation* animation)
{
	bov_points_delete(animation->particles_flow);
    bov_points_delete(animation->particles_pool);
	if(animation->grid != NULL)
		bov_points_delete(animation->grid);
	bov_window_delete(animation->window);
	free(animation);
}
void display_particles_boundary(Particle** particles, Animation* animation,bool end, int iter, double bounds[6]){
	int N_pool = animation->N_pool;
    int N_flow = animation->N_flow;
	GLfloat(*data_pool)[8] = malloc(sizeof(data_pool[0])*N_pool);
    GLfloat(*data_flow)[8] = malloc(sizeof(data_flow[0])*N_flow);
	fillData(data_pool, data_flow, particles, N_pool, N_flow);
// 	colours_neighbors(data, particles, N / 2);
	animation->particles_pool = bov_particles_update(animation->particles_pool,data_pool,N_pool);
    animation->particles_flow = bov_particles_update(animation->particles_flow,data_flow,N_flow);
	free(data_pool);
    free(data_flow);

	char screenshot_name[64] = "dambreak_";
	char int_string[32];
	sprintf(int_string, "%d", iter);
	strcat(screenshot_name, int_string);
 	strcat(screenshot_name, ".png");
	//DRAW BOUNDARIES
	double xl1 = bounds[0];
	double xl2 = bounds[1];
	double xr = bounds[2];
	double yb1 = bounds[3];
    double yb2 = bounds[4];
    double yt = bounds[5];
	GLfloat coord_bound[][2] = {
        { xl1, yb1},
        { xl2, yb1},
        { xl2, yb2},
        { xr, yb2},
        { xr,  yt},
        { xl1,  yt},
        { xl1,  yb1},
	};
	bov_points_t* point_bound = bov_points_new(coord_bound, 6, GL_STATIC_DRAW);
	bov_points_param_t BoundParams = {
		.fillColor = { 0.0, 1.0, 0.0, 1.0 },
		.scale = { 1.0, 1.0 },
		.width = 0.001
	};
	bov_points_set_param(point_bound, BoundParams);
	
	bov_window_t* window = animation->window;
    double tbegin = bov_window_get_time(window);
	if (!end){
		while (bov_window_get_time(window) - tbegin < animation->timeout) {
			if(animation->grid != NULL)
                bov_lines_draw(window,animation->grid,0, BOV_TILL_END); //To show the grid
				bov_line_loop_draw(window, point_bound, 0, BOV_TILL_END);
			bov_particles_draw(window, animation->particles_flow, 0, BOV_TILL_END);
            bov_particles_draw(window, animation->particles_pool, 0, BOV_TILL_END);
            //display_neighbours(animation->window, animation, particles, N_flow);
            //printf("no problem\n");
			if (iter%1000 == 0) bov_window_screenshot(window, screenshot_name);
			bov_window_update(window);
		}
	}
	else {
		// we want to keep the window open with everything displayed
		while (!bov_window_should_close(window)) {
			if (animation->grid != NULL)
				bov_lines_draw(window, animation->grid, 0, BOV_TILL_END);
			bov_particles_draw(window, animation->particles_flow, 0, BOV_TILL_END);
            bov_particles_draw(window, animation->particles_pool, 0, BOV_TILL_END);
			bov_window_screenshot(window, screenshot_name);
			bov_window_update_and_wait_events(window);
		}
	}
    bov_points_delete(point_bound);
}

/*
void display_particles(Particle** particles, Animation* animation,bool end, int iter)
{
	int N = animation->N;
	GLfloat(*data)[8] = malloc(sizeof(data[0])*N);
	fillData(data, particles, N);
// 	colours_neighbors(data, particles, N / 2);
	animation->particles = bov_particles_update(animation->particles,data,N);
	free(data);

	char screenshot_name[64] = "circle_";
	char int_string[32];
	sprintf(int_string, "%d", iter);
	strcat(screenshot_name, int_string);
// 	strcat(screenshot_name, ".png");
	bov_window_t* window = animation->window;
	double tbegin = bov_window_get_time(window);
	if (!end){
		while (bov_window_get_time(window) - tbegin < animation->timeout) {
			if(animation->grid != NULL)
				bov_lines_draw(window,animation->grid,0, BOV_TILL_END); //To show the grid
			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
            display_neighbours(animation->window, animation, particles, N);
			//if (iter%10 == 0) bov_window_screenshot(window, screenshot_name);
			bov_window_update(window);
		}
	}
	else {
		// we want to keep the window open with everything displayed
		while (!bov_window_should_close(window)) {
			if (animation->grid != NULL)
				bov_lines_draw(window, animation->grid, 0, BOV_TILL_END);
			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
			bov_window_screenshot(window, screenshot_name);
			bov_window_update_and_wait_events(window);
		}
	}
}
 */


void display_neighbours(bov_window_t* window, Animation* animation, Particle** particles, int N)
{
    int i;
    int j;
    Particle* p;
    int num_neighbours;
    for (i=2800; i<=2800; i++) {
        j = 0;
        num_neighbours = 0;
        p = particles[i];
        ListNode* node = p->neighborhood->head;
        while(node != NULL){
            num_neighbours += 1;
            node = node->next;
        }
        GLfloat(*neighbours)[2] = malloc(sizeof(neighbours[0]) * 2 * num_neighbours);
        node = p->neighborhood->head;
        while (node != NULL) {
            Particle* q = (Particle*)node->v;
            neighbours[j][0] = q->pos->x;
            neighbours[j][1] = q->pos->y;
            j+=1;
            node = node->next;
        }
        Node_free(node);
        bov_points_t *points = bov_points_new(neighbours, num_neighbours , GL_STATIC_DRAW);
        GLfloat colors[] = { 0.0, 1.0, 0.0, 1.0 };
        bov_points_set_color(points, colors);
        bov_points_set_width(points, 0.001);
        bov_points_set_outline_width(points, 0.0001);
        bov_points_draw(window, points, 0, BOV_TILL_END);
        bov_points_delete(points);
        free(neighbours);
    }
}






// void display_particles(Particle** particles, Animation* animation,bool end)
// {
// 	int N = animation->N;
// 	GLfloat(*data)[8] = malloc(sizeof(data[0])*N);
// 	fillData(data, particles, N);
// // 	colours_neighbors(data, particles, N / 2);
// 	animation->particles = bov_particles_update(animation->particles,data,N);
// 	free(data);
//
// 	bov_window_t* window = animation->window;
// 	double tbegin = bov_window_get_time(window);
// 	if (!end){
// 		while (bov_window_get_time(window) - tbegin < animation->timeout) {
// 			if(animation->grid != NULL)
// 				bov_lines_draw(window,animation->grid,0, BOV_TILL_END);
// 			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
// 			bov_window_update(window);
// 		}
// 	}
// 	else {
// 		// we want to keep the window open with everything displayed
// 		while (!bov_window_should_close(window)) {
// 			if (animation->grid != NULL)
// 				bov_lines_draw(window, animation->grid, 0, BOV_TILL_END);
// 			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
// 			bov_window_update_and_wait_events(window);
// 		}
// 	}
// }

// Fills color with the color to use for particle p
void colormap_cell(Particle* p, float color[3]) {
	if (p->cell == NULL) {
		color[0] = 0;color[1] = 0;color[2] = 0;
	}
	else {
		if (p->cell->i % 2 == 0) {
			color[0] = 20;
			if (p->cell->j % 2 == 0)
				color[1] = 20;
			else
				color[1] = 0;
		}
		else {
			color[0] = 0;
			if (p->cell->j % 2 == 0)
				color[1] = 20;
			else
				color[1] = 0;
		}
		color[2] = 0;
	}
}

void colormap_uni_color(float color[3])
{
	color[0] = 0;color[1] = 10;color[2] = 20;

}

void colormap_uni_color_2(float color[3]) {
	color[0] = 0;color[1] = 0;color[2] = 20;
}

void colormap_Cs(Particle *p, float color[3]) {
	color[0] = 20*squared(p->Cs);
	color[1] = 0;
	color[2] = 20*squared(1.0-p->Cs);
}

void colormap_fs(Particle *p, float color[3], double max_norm) {
	xy* fs = xy_new(-p->param->sigma * p->normal->x * p->kappa, -p->param->sigma * p->normal->y * p->kappa);
// 	double fs_norm = 1.0;//norm(fs);
	color[0] = 10*squared(fs->x/max_norm) + 10*squared(fs->y/max_norm);
	color[1] = 0.0;
	color[2] = 0.0;//20*squared(fs->y/max_norm);
    free(fs);
}

void colours_neighbors(GLfloat(*data)[8], Particle** particles, int index)
{
	Particle* p = particles[index];
	ListNode *node = p->neighborhood->head;
	data[index][4] = 20;data[index][5] = 20;data[index][6] = 20;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		int i = q->index;
		data[i][4] = 0;data[i][5] = 20;data[i][6] = 20;
		node = node->next;
	}
    Node_free(node);
}
