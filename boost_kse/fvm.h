#ifndef __FVM_H__
#define __FVM_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kse.h"

extern double gamma_g;

void execute(State &st, prev_layer *pl);

#endif
