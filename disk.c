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
    if (fabs(VISCOSITYRATIO-1.0)<0.1) return viscosity;
		if (r < rmin) viscosity *=  VISCOSITYRATIO;
		if ((r >= rmin) && (r <= rmax)) {
    viscosity *= 0.5*(1.0-VISCOSITYRATIO)*sin((r-r_tran)*M_PI/(2.0*dr_tran))+0.5*(1.0+VISCOSITYRATIO);
		}

		double c0,c1,c2,c3,rlog,alplog,rtran1,rtran2,rtran0;
		rlog=log10(r);
		if(MDOT_INT==0) {
		rtran0=0.064;rtran1=0.205;rtran2=0.73;
		if(r<rtran0)	alplog=log10(0.09848);
		else if(r>=rtran0 && r < rtran1){
			c0=-6.7811175; c1=-12.02901; c2=-8.2440562; c3=-1.8594328;
			alplog=c0+c1*rlog+c2*rlog*rlog+c3*rlog*rlog*rlog;
		}
		else if(r>=rtran1 && r<rtran2){
			c0=-4.5536009; c1=-4.1222938; c2=-0.17861849;
			alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}else if(MDOT_INT==1){
		rtran0=0.014;rtran1=0.049;rtran2=0.232;
		if(r<rtran0)  alplog=log10(0.09848);
		else if(r>=rtran0 && r < rtran1){
		c0=-8.5526396; c1=-9.0739395; c2=-3.1434803; c3=-0.22451420;
		alplog=c0+c1*rlog+c2*rlog*rlog+c3*rlog*rlog*rlog;
													    }
		else if(r>=rtran1 && r<rtran2){
		c0=-6.6333076; c1=-4.3837696; c2=-0.37942259;
		alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}else if(MDOT_INT==2){
		rtran2=0.066;
		if(r<rtran2){
		c0=-11.604954; c1= -9.0560947;c2= -2.2069353;
		alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}
		if (SINEALPHA==1) return viscosity;
		else return pow(10,alplog);
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
if ( (mdot<0e-10 && alpha>20e-4) || r/LUNIT>10.0 || temper_passive > temper_active) return temper_passive;
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
        double dr,i_r,r1,r2;
        int i;
        i_r=ring_num*log(r/r_min)/log(R_OUT/r_min);
        i=floor(i_r);
        //if(!ITER) printf("k_P func, r=%f\ti_r=%f\n",r,i_r);
        dr=0.002;
        r1=r_min*exp(i*1.0/ring_num*log(R_OUT/r_min));
        r2=r_min*exp((i+1)*1.0/ring_num*log(R_OUT/r_min));
        return -1.0*(log(pressure(r1))-log(pressure(r2)))/(log(r1)-log(r2));
        //return -1.0*(log(pressure(r))-log(pressure(r+dr)))/(log(r)-log(r+dr));
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
