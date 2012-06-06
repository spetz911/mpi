#ifndef __FVM_H__
#define __FVM_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kse.h"


void init_state(state*, double, double, double);
void execute(state*, prev_layer*, double);

#endif
