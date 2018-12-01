#include <stdio.h>
#include "ex_func.h"
#include "global_ex.h"

double vr_p[2]={0.0};
void check_disk(double r){
	FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9;
	int i,j;
	ITER=0;
	opa=10.0;
	printf("Properties at %1.2fAU:mdot=%e alpha=%f\n",r,mdot,alpha);
	printf("Surf_dens=%g\nTemperature=%f\nMidplane-density=%g\nSound_speed=%g\n\
H/R=%g\nMean free path=%g\nv_K=%g\nvt_gas=%g\nvr_peb=%g\nzeta=%g\n\
Omega_K=%g\nvr_gas=%g\nmdot=%g\n",\
Sigma(r),temperature(r),density(r)\
,sound_sp(r),height(r)/r/LUNIT,mean_path(r),v_K(r),\
vt_gas(r),vr_estimate(1.0,1.0,vr_p),yeta(r),w_K(r),\
vr_gas(r),vr_gas(r)*2*M_PI*r*LUNIT*Sigma(r)*TUNIT/MUNIT);

	fp1=fopen("opacity.txt","w");
	fp2=fopen("sigma.txt","w");
	fp3=fopen("temperature.txt","w");
	fp4=fopen("dust_mass.txt","w");
	fp5=fopen("height_peb.txt","w");
	fp6=fopen("vr_peb.txt","w");
	fp7=fopen("alpha.txt","w");
	fp8=fopen("vr_gas.txt","w");
	fp9=fopen("dpdr.txt","w");
	for(i=0;i<ring_num;i++){
	fprintf(fp1,"%e\t%e\n",dust_budget[i].rad,func_line1(dust_budget[i].rad,p_opa_line));
	fprintf(fp2,"%e\t%e\n",dust_budget[i].rad,Sigma(dust_budget[i].rad));
	fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,temperature(dust_budget[i].rad));
	fprintf(fp4,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].mass_out);
	fprintf(fp7,"%e\t%e\n",dust_budget[i].rad,alpha_func(dust_budget[i].rad));
	fprintf(fp8,"%e\t%e\n",dust_budget[i].rad,vr_gas(dust_budget[i].rad));
	fprintf(fp9,"%e\t%e\n",dust_budget[i].rad,k_P_func(dust_budget[i].rad));

	for(j=0;j<peb_size_num;j++){
	fprintf(fp5,"%e\t",peb_map[i].hei[j]);
	fprintf(fp6,"%e\t",peb_map[i].vr[j]);
	}
	fprintf(fp5,"\n");
	fprintf(fp6,"\n");
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);
	fclose(fp7);
	fclose(fp8);
	fclose(fp9);
}

