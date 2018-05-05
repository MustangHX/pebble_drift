//
//  main.c
//  pebble_size_drift
//
//  Created by Xiao Hu on 5/16/15.
//
//
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
double pp_vr_tau0[2]={0.0};
double pp_vr_tau2[2]={0.0};
int main(int argc, char *argv[])
//	int argc;
//	char *argv[];
{
	alpha=alpha_init;
	mdot=mdot_init;
	time_t time0,time1;
	int i,j,n,num_step=0,tot_num_step=(int)(time_yr*1.0/init_step),check,NbRestart;
	double AREA,a_pb1,a_max,vol_plus,tau,vr0,mass_flow_inner;
	double coag_eff=1.0,tot_mass=0.0,out_source=0.0,a_p,r0,dt=init_step,time_sum=0.0,dt2,tot_mass_dust,t_single=0.0;
	double vr1,vr2,tau1,tau2;
	FILE *fp,*fp2,*fp3,*fp4;
	char outname[256], outname2[256];
	int Restarting = 0,grow_method=3;
	printf("%s\n",argv[0]);
	for(i=0; i< argc; i++){
		printf("%s",argv[i]);
		if (strchr (argv[i], 's')) {
		Restarting = 1;
		i++;
		NbRestart = atoi(argv[i]);
		}
		if (NbRestart < 0) printf("Incorrect restart number\n");
	}
	n=0;
   // group();
	printf("PPPPeEEEBBBB%.16f\n",density(1.0));
	//drift(1.0,drag_group(1.0,0.01));
        printf("PPPPeEEEBBBB%.12f\n",drag_group(1.0,0.01));
	Init2();
	printf("SIZEOF PEB=%lu\t DUST=%lu\n",sizeof(peb_map[0]),sizeof(dust_budget[0]));
	check_disk(1.0);
	drift_vr_test(1.0,152.0,pp_vr_tau2);
        fp=fopen("vr_check.txt","w");
	printf("===============================\n");
	double re1;
	for(i=1;i<3000;i++){
	a_pb1=i*0.1;
	r0=1.0;
	//a_pb1=9.0;
        vr1=vr_estimate(r0,a_pb1,pp_vr_tau0);
        tau1=pp_vr_tau0[1];
        vr2=drift_vr(r0,a_pb1,pp_vr_tau2);
        tau2=pp_vr_tau2[1];
	re1=2.0*a_pb1*sqrt(vr1*vr1+0.25*tau1*tau1*vr1*vr1)/viscosity(r0);
	//printf("SIZE=%f\n",a_pb1);
        fprintf(fp,"%f %f\t est= %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",2.25*mean_path(r0),a_pb1,vr1,2.0*a_pb1*sqrt(vr1*vr1+0.25*tau1*tau1*vr1*vr1)/viscosity(r0),\
                tau1,vr2,2.0*a_pb1*sqrt(vr2*vr2+0.25*tau2*tau2*vr2*vr2)/viscosity(r0),tau2,v_r2(r0,a_pb1)[0],v_r2(r0,a_pb1)[1]);
       	//printf("%f\t%f\t%f\n",vr1,v_r1_test(1.0,a_pb1),re1);
	}
	//printf("%f\t%f\n",vr_estimate(1.0,a_pb1,pp_vr_tau0),v_r1_test(1.0,a_pb1));
        fclose(fp);
	printf("vr_check finished\n");
	//return 0;
	if(Restarting == 1){
		Restart(NbRestart);
		time_sum+=NbRestart;
	}
        fp=fopen("size_chart.txt","w");
        for(i=0;i<peb_size_num;i++){
                fprintf(fp,"%f\n",peb_map[0].size_med[i]);
	        }
        fclose(fp);
	        fp=fopen("rad_chart.txt","w");
        for(i=0;i<ring_num;i++){
                fprintf(fp,"%f\n",peb_map[i].rad);
	        }
				fprintf(fp,"%f\n",R_OUT);
        fclose(fp);
	printf("disk_grid_num finished\n");
	fp=fopen("tau_chart.txt","w");
	for(i=0;i<ring_num;i++){
	for(j=0;j<peb_size_num;j++){
		fprintf(fp,"%e\t",peb_map[i].tau_fric[j]);
	}
	fprintf(fp,"\n");
	}
	fclose(fp);
	out_source=0.0;
	for(j=0;j<peb_size_num;j++){
		i=ring_num-1;
		AREA=peb_map[i].AREA;
		if(peb_map[i].size_med[j]>size_min_inj && peb_map[i].size_med[j]<peb_size_lim){
		out_source+=AREA*p_size_func(peb_map[i].size_med[j]);
		}
	}
	printf("inject_normalize finished\n");
	fp=fopen("1mm.txt","w");
	num_step=0;
	a_p=0.1;
	r0=30.00;
        fprintf(fp,"%e\%e\n",r0,a_p);

	while(0 && r0>0.1){
	vr0=drift_vr(r0,a_p,pp_vr_tau0);
	tau=pp_vr_tau0[1];
	r0-=vr0*dt*TUNIT/LUNIT;
	t_single+=dt;
	a_max=2.25*mean_path(r0);
	vol_plus=0.0*M_PI*a_p*a_p*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT;
	if(a_p > a_max) {check=1; vol_plus=0.0;}
	else check=0;
	a_p=pow(((vol_plus*coag_eff*dust_gas*density(r0)/rho_peb+4.0/3.0*M_PI*a_p*a_p*a_p)*3.0/4.0/M_PI),0.33333333333333333);
	if(check==0 && a_p>=a_max)       a_p=a_max;
	fprintf(fp,"%e\t%e\t%e\n",r0,a_p,t_single);
	}
	fclose(fp);
	printf("single_integration finished\n");
	fp=fopen("gas_vr.txt","w");
	for(i=0;i<ring_num;i++){
	fprintf(fp,"%f\t%g\t%g\n",dust_budget[i].rad_med,vr_gas(dust_budget[i].rad_med),vr_gas(dust_budget[i].rad_med)*TUNIT/LUNIT);
	}
	fclose(fp);
	fp=fopen("gas_dust_density.txt","w");
        for(i=0;i<ring_num;i++){
        fprintf(fp,"%f\t%g\t%g\n",dust_budget[i].rad_med,Sigma(dust_budget[i].rad_med),dust_budget[i].surf_dens);
        }
        fclose(fp);
	printf("gas_dust_dens finished\n");
	stokes_size();
        printf("stokes_size finished\n");
	if(!(mdot<2e-10 && alpha>8e-4)) tau_unity();
        printf("tau_unity finished\n");
	num_step=0;
	if(Restarting == 1){
		num_step=NbRestart;
	}

	time0=time(NULL);
	//start time sequence
	//disk_evolve();
	printf("main integration started\n");
	while (time_sum<tot_num_step*1.0)
	{
	if(num_step==0){
        sprintf(outname,"out_sigma%d.txt",num_step);
        fp=fopen(outname,"w");
        for(i=0;i<ring_num;i++){
        for(j=0;j<peb_size_num;j++){
        fprintf(fp,"%e\t",peb_map[i].surf_dens[j]);
	        }
        fprintf(fp,"\n");
	}
        fclose(fp);
	sprintf(outname2,"dust_sigma%d.txt",(int)time_sum);
	fp3=fopen(outname2,"w");
	for(i=0;i<ring_num;i++){
		tot_mass+=dust_budget[i].mass_out;
		fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].surf_dens);
	}
	fclose(fp3);
	}
	if(grow_method==3){
		dt2=dt;
	if (1 && ((int)time_sum)%1000==0){
		alpha=alpha_init;
		mdot=mdot_init*exp(-1.0*time_sum/1e6*log(1e-8/1e-9));
		//check_disk(1.0);
		//disk_evolve();
		printf("disk evolved\n");
	}
	dt=grow_3b_ada_fix(dt2,time_sum);
        //      dt=grow_3b2_test(dt2);

		//printf("main dt=%f\tdt2=%f\n",dt,dt2);
		//dt=1.0;
		if(dt<0.0) {
		printf("Actual time step count:%d\t dt=%f\n",num_step,dt);
		return 0;
		}
	}
	//dust_evolve(dt);
	//disk_evolve();
	if(COAG_SW>0) {coagulation(dt,time_sum);}
	mass_flow_inner=0.0;
	for(i=ring_num-1;i>-1;i--){
        for(j=0;j<peb_size_num;j++){
                if(i==0) mass_flow_inner+=peb_map[i].surf_dens[j]*peb_map[i].AREA;

		peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
		peb_map[i].mass_in[j]=0.0;
		//AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
		AREA=peb_map[i].AREA;

		if(i==ring_num-1 && 1 && 1) {
			//peb_map[i].mass_out[j]=0.2*AREA*dust_budget[i].surf_dens*exp(-1.0*peb_map[i].size[j]/0.1)*exp(0.0*num_step/100);
		if(peb_map[i].size_med[j]>size_min_inj && peb_map[i].size_med[j]<peb_size_lim){
			peb_map[i].mass_out[j]+=1.0*dust_gas*peb_dust_inj*AREA*MUNIT*mdot*dt*p_size_func(peb_map[i].size_med[j])/out_source;
			}
		}
		else if(i<ring_num-1 && j<10 && 0){
                        peb_map[i].mass_out[j]=0.1*AREA*dust_budget[i].surf_dens*p_size_func(peb_map[i].size_med[j])*exp(0.0*num_step/100);
		}

		else peb_map[i].surf_dens[j]+=peb_low_lim/10;
		peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
		peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
		//printf("peb_rho=%e\t%e\t",peb_map[i].rho[j],peb_map[i].surf_dens[j]);
	}
	}
	//frag();
	num_step++;
	time_sum+=dt;
	if(time_sum-floor(time_sum)<0.00001 && ((int)(time_sum))%((int)(outp_step))==0){
		tot_mass=0.0;
		tot_mass_dust=0.0;
        sprintf(outname,"out_sigma%d.txt",(int)time_sum);
        fp=fopen(outname,"w");
	sprintf(outname2,"dust_sigma%d.txt",(int)time_sum);
	fp3=fopen(outname2,"w");
	fp4=fopen("mass_flow_inner_ring.txt","a+");
	fp2=fopen("mass_check.txt","a+");
        for(i=0;i<ring_num;i++){
        for(j=0;j<peb_size_num;j++){
		if(peb_map[i].surf_dens[j] < peb_low_lim) peb_map[i].surf_dens[j]=peb_low_lim;
        fprintf(fp,"%e\t",peb_map[i].surf_dens[j]);
		tot_mass+=peb_map[i].mass_out[j];
	}
        fprintf(fp,"\n");
	}
	for(i=0;i<ring_num;i++){
		tot_mass_dust+=dust_budget[i].mass_out;
	fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].surf_dens);
	}
	fprintf(fp2,"%2.20g\t%2.20g\t%2.20g\n",tot_mass,tot_mass_dust,tot_mass+tot_mass_dust);
	fprintf(fp4,"%f\t%2.20g\n",time_sum,mass_flow_inner);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	printf("%f finished\r",time_sum/(tot_num_step*1.0));
	printf("Actual time step count:%d\t dt=%f\t time=%f\n",num_step,dt,time_sum);
	}

	}
	time1=time(NULL);

	printf("Actual time step count:%d dt=%f time_used=%f min\n",num_step,dt,(time1-time0)/60.0);

return 0;
}
