#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#include "global_ex.h"

double v_K(double r){
    return 29.8*100000.0*sqrt(m_star)/sqrt(r);
}

double w_K(double r){
    return v_K(r)/(r*LUNIT);
}

double alpha_func(double r){
	  double rmin = RTRAN-DRTRAN;
	  double rmax = RTRAN+DRTRAN;
	  double dr_tran = DRTRAN;
	  double r_tran = RTRAN;
		double viscosity=alpha_init;
		if (r < rmin) viscosity *=  VISCOSITYRATIO;
		if ((r >= rmin) && (r <= rmax)) {
    viscosity *= 0.5*(1.0-VISCOSITYRATIO)*sin((r-r_tran)*M_PI/(2.0*dr_tran))+0.5*(1.0+VISCOSITYRATIO);
		}
		return viscosity;
}

double temperature (double r){
alpha = alpha_func(r);
r=r*LUNIT;
double temper_active,temper_passive;
if (!ITER) {
	opa=func_line1(r/LUNIT,p_opa_line);
//	printf("rad=%e\tOPA=%e\n",r,opa);
}
 temper_active=pow(3.0,0.2)*pow(2.0,-1.4)*pow(M_PI,-0.4)\
        *pow(mu*m_p/gamma0/k_B,0.2)*pow(opa/sig_sb,0.2)\
	*pow(alpha,-0.2)*pow(G*m_star*MUNIT,0.3)\
	*pow((1-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.4)*pow(r,-0.9);

temper_passive=temp0*pow(r/LUNIT,-3.0/7.0);
if ( (mdot<2e-10 && alpha>8e-4) || r/LUNIT>10.0 || temper_passive > temper_active) return temper_passive;
else return temper_active;
}	

double sound_sp(double r) {//return c_s in cgs
if (!ITER) opa=func_line1(r,p_opa_line);
return sqrt(gamma0*k_B*temperature(r)/mu/m_p);
}

double height(double r) {// return scale height in cgs
if (!ITER) opa=func_line1(r,p_opa_line);
return sound_sp(r)/w_K(r);
}

double Sigma (double r) {//return surface density in cgs
alpha = alpha_func(r);
if (!ITER) opa=func_line1(r,p_opa_line);
return mdot*MUNIT/TUNIT/3.0/M_PI/(alpha*sound_sp(r)*height(r));
}

double density(double r) {
if (!ITER) opa=func_line1(r,p_opa_line);
return Sigma(r)/height(r)/sqrt(2*M_PI);
}

double mean_path(double r){

if (!ITER) opa=func_line1(r,p_opa_line);
return mu*m_p/density(r)/2e-15;
}

double viscosity( double r){//http://www.ifu.ethz.ch/IE/education/AQAM/GasKinetics
	    //return 0.4991*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r)*100000.0;
if (!ITER) opa=func_line1(r,p_opa_line);
	return 0.5*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r);
}

double pressure(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
	return temperature(r)*k_B*density(r)/mu/m_p;
}

double k_P_func(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
	double dr;
	dr=0.002;
	return -1.0*(log(pressure(r))-log(pressure(r+dr)))/(log(r)-log(r+dr));
}	

double yeta(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
    return k_P_func(r)*(sound_sp(r)/v_K(r))*(sound_sp(r)/v_K(r));
}

double vt_gas(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
        return v_K(r)*sqrt(1-yeta(r));
}


double vr_gas (double r){
alpha = alpha_func(r);
if (!ITER) opa=func_line1(r,p_opa_line);
return 3.0*alpha*sound_sp(r)*height(r)/2.0/r/LUNIT;
}
