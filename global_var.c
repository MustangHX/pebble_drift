#include "global_var.h"
PEBBLE peb_group[peb_num];
PEBBLE_MAP peb_map[ring_num];
DUST_MAP dust_budget[ring_num];
double dt_ring[ring_num];
double dtco_ring[ring_num];
double sigma_pre[ring_num];
double sigma_aft[ring_num];
double alpha;
double mdot;
double opa;
int ITER;
SPLINE opa_line;

pSPLINE p_opa_line = &opa_line;
