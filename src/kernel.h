#ifndef KERNEL_H
#define KERNEL_H
#include "utils.h"

typedef enum Kernel Kernel;

enum Kernel {Cubic};

xy* grad_kernel(xy* p1, xy* p2, double kh, Kernel kernel);
double eval_kernel(xy *p1, xy *p2, double kh, Kernel kernel);
double deriv2_Cubic_kernel(xy* p1, xy* p2, double kh, Kernel kernel);

#endif
