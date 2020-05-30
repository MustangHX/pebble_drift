#include "global_ex.h"
#include "global_var.h"
#include "ex_func.h"
#include<math.h>
#include<stdio.h>
double pp_vr_tau[2]={0.0};


double grow_3b_ada_fix(double dt0, double tot_time){ //adaptive timestep with variable radial resolution
        int i,j,i_new,j_new;
        double a_pb1,a_pb11,a_pb2,a_pb22,a_pb3,vr0,vr1,vr2,vt0,AREA,dr,a_max,dt1,sub_time1,mass_gain,ring_mass_gain=0.0,ring_sigma;
        double tau,vol_plus,frac,frac_s,coag_eff,ring_mass_before,ring_mass_after,old_sigma,ratio_sigma=1.0,rho_eff,h1,h2;
        dr=size_ring;
        coag_eff=0.0;


if(1 && ((int)tot_time)%100==0 || tot_time<10000.0){
for(i=0;i<ring_num;i++){
      dt_ring[i]=dt0;
        }
}

for(i=ring_num-1;i>=0;i--){
for(j=0;j<peb_size_num;j++){
	peb_map[i].fluxL[j]=0.0;
	peb_map[i].fluxR[j]=peb_map[i].fluxRb[j];
	peb_map[i].fluxRb[j]=0.0;
}
}
for(i=ring_num-1;i>=0;i--){
AREA=peb_map[i].AREA;
dr=peb_map[i].dr;
a_max=2.25*mean_path(peb_map[i].rad+dr/2.0);
dt1=2.0*dt0;

if(1 && ((int)tot_time)%100!=0 && tot_time>50000.0 ){
	         dt1=dt_ring[i];
}
else {
	dt1=2.0*dt0;
do{
dt1=dt1/2.0;
mass_gain=0.0;
frac=0.0;
for(j=0;j<peb_size_num;j++){
	if(1 && peb_map[i].surf_dens[j]<peb_low_lim*1e10) continue;
        a_pb1=peb_map[i].size[j];
        //vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb1,pp_vr_tau);
	vr0=peb_map[i].vr_med_r[j];
	vt0=peb_map[i].vt_med_r[j];
	h1=dust_budget[i].hei;
	h2=peb_map[i].hei[j];
        if(fabs(vr0*dt1*TUNIT/LUNIT/dr)>fabs(frac)) frac=vr0*dt1*TUNIT/LUNIT/dr;
	rho_eff=dust_budget[i].surf_dens/sqrt(2*M_PI*(h1*h1+h2*h2));
	tau=pp_vr_tau[1];
        if(a_pb1>a_max ||(tot_time>1e4 && i==0)) vol_plus=0.0;
        else{
        vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+vt0*vt0)*dt1*TUNIT;
        }
        a_pb11=pow(((vol_plus*coag_eff*rho_eff/rho_peb+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),1.0/3.0);
        if(a_pb1<a_max && a_pb11>a_max) a_pb11=a_max;
	if(a_pb11/a_pb1-1.0>mass_gain) mass_gain=a_pb11/a_pb1-1.0;
}

//printf("mass_gain=%f\t %f\t %f\t ring_num=%d\n",mass_gain,vol_plus,dt1,i);
}while(mass_gain>0.003 || fabs(frac)>0.09);
dt_ring[i]=dt1;
}

//dt1=dt0;
//if(i<i_lim1) dt1=dt0/10.0;
if(0 && dt1<dt0) printf("GROWTH:RING_NUM=%d\tnew_dt=%f\n",i,dt1);
sub_time1=0.0;
//if(i==0) dt1=dt1/40.0;
while(sub_time1<dt0-dt1/10.0){
//printf("SUB_TIME=%f\n",sub_time1);
ring_mass_before=0.0;
ring_mass_after=0.0;
ring_mass_gain=0.0;
ring_sigma=0.0;
for(j=0;j<peb_size_num;j++){
        ring_mass_before+=peb_map[i].mass_out[j];
	ring_sigma+=peb_map[i].surf_dens[j];
}
for(j=0;j<peb_size_num;j++){
        a_pb1=peb_map[i].size[j];
        //vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb1,pp_vr_tau);
  	vr0=peb_map[i].vr_med_r[j];
        vt0=peb_map[i].vt_med_r[j];

	h1=dust_budget[i].hei;
	h2=peb_map[i].hei[j];
	rho_eff=dust_budget[i].surf_dens/sqrt(2*M_PI*(h1*h1+h2*h2));
	//tau=pp_vr_tau[1];
	
        if(a_pb1>a_max ||dust_budget[i].surf_dens< 1e-6 || fabs(ring_sigma)/dust_budget[i].surf_dens>1e5) vol_plus=0.0;
        //if(a_pb1>a_max ||dust_budget[i].surf_dens< 1e-6 || i==0) vol_plus=0.0;
        else{
        //vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt1*TUNIT;
	//printf("vol1=%g\t",vol_plus);
        vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+vt0*vt0)*dt1*TUNIT;
	//printf("vol2=%g\n",vol_plus);
        }
        a_pb11=pow(((vol_plus*coag_eff*rho_eff/rho_peb+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),1.0/3.0);
        if(a_pb1<a_max && a_pb11>a_max) a_pb11=a_max;
        a_pb2=peb_map[i].size[j+1];
       // vr0=vr_estimate(peb_map[i].rad+dr/2.0,a_pb2,pp_vr_tau);
	vr0=peb_map[i].vr_med_r[j+1];
        vt0=peb_map[i].vt_med_r[j+1];

        if(a_pb2>a_max || dust_budget[i].surf_dens< 1e-6 || fabs(ring_sigma)/dust_budget[i].surf_dens>1e5) vol_plus=0.0;
        //if(a_pb2>a_max || dust_budget[i].surf_dens< 1e-6 || i==0) vol_plus=0.0;
        else{
        //    vol_plus=1.0*M_PI*a_pb2*a_pb2*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt1*TUNIT;
	vol_plus=1.0*M_PI*a_pb2*a_pb2*sqrt(vr0*vr0+vt0*vt0)*dt1*TUNIT;
        
	    }
        a_pb22=pow(((vol_plus*coag_eff*rho_eff/rho_peb+4.0/3.0*M_PI*a_pb2*a_pb2*a_pb2)*3.0/4.0/M_PI),1.0/3.0);
        if(a_pb2<a_max && a_pb22>a_max) a_pb22=a_max;
	
        a_pb3=(a_pb2-a_pb1)/(a_pb22-a_pb11)*(a_pb2-a_pb11)+a_pb1;
        frac_s=(a_pb2-a_pb3)/(a_pb2-a_pb1);
//      frac_s=(a_pb22-a_pb2)/(a_pb2-a_pb1);
        if(a_pb2>a_max) frac_s=0.0;
	//ring_mass_gain+=frac_s*pow(a_pb2/a_pb1,3.0)*peb_map[i].mass_out[j];
	ring_mass_gain+=(pow(a_pb2/a_pb3,3.0)-1.0)*peb_map[i].mass_out[j];

        if(frac_s>0.5 && 0) printf("%f\t %d\t%d\t sizeTOOLARGE_middle\n",frac_s,i,j);
	if(a_pb1>a_max && frac_s>0.0) printf("WTF??? a_max=%f\ta_pb1=%f\tfrac_s=%f\n",a_max,a_pb1,frac_s);
        //vr1=vr_estimate(peb_map[i].rad+dr,(a_pb1+a_pb2)/2.0,pp_vr_tau);
        //vr2=vr_estimate(peb_map[i].rad,(a_pb1+a_pb2)/2.0,pp_vr_tau);
  if(i<ring_num-1) vr1=peb_map[i+1].vr_med_s[j];
	else  vr1=vr_estimate(peb_map[i].rad+dr,(a_pb1+a_pb2)/2.0,pp_vr_tau);
				vr2=peb_map[i].vr_med_s[j];
			  frac=1.0-(dr-vr2*dt1*TUNIT/LUNIT)/(dr+vr1*dt1*TUNIT/LUNIT-vr2*dt1*TUNIT/LUNIT);
				frac+=peb_map[i].vr_drag[j]*dt1*TUNIT/LUNIT/dr;

        j_new=j+1;
	if(fabs(frac)>0.2 && 0) printf("OMG moving too fast%d\t%d\t%f\n",i,j,frac);
        if(a_pb2>a_max || j_new>peb_size_num-1 ) {j_new=j; frac_s=0.0;}
        i_new=i-1;
//        if(i_new<0) i_new=0;
	if(1 && peb_map[i].surf_dens[j]<peb_low_lim*1e10){
	frac=0.0;
	frac_s=0.0;
	continue;
	}
/*	peb_map[i_new].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*(1.0-frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i_new].mass_in[j]+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j]+=pow(a_pb2/a_pb3,3)*(1.0-frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
*/
	if(1 && i_new<0 && frac>0.0 ){
	i_new=i;
	peb_map[i_new].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*fabs(frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i_new].mass_in[j]+=pow(a_pb2/a_pb3,3)*fabs(frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	}
	if(1 && i==ring_num-1 && frac<0.0 ){
	peb_map[i].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*fabs(frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j]+=pow(a_pb2/a_pb3,3)*fabs(frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	}
	if(frac>=0.0){
	peb_map[i].fluxL[j_new]+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	peb_map[i].fluxL[j]+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	}
	else{
	peb_map[i].fluxRb[j_new]+=pow(a_pb22/a_pb2,3)*fabs(frac)*frac_s*peb_map[i].mass_out[j];
	peb_map[i].fluxRb[j]+=pow(a_pb2/a_pb3,3)*fabs(frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
	}

	peb_map[i].mass_in[j_new]+=pow(a_pb22/a_pb2,3)*(1.0-fabs(frac))*frac_s*peb_map[i].mass_out[j];
	peb_map[i].mass_in[j]+=pow(a_pb2/a_pb3,3)*(1.0-fabs(frac))*(1.0-frac_s)*peb_map[i].mass_out[j];
	peb_map[i].mass_out[j]-=peb_map[i].mass_out[j];
	
/*
	ring_mass_after+=pow(a_pb22/a_pb2,3)*frac*frac_s*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb22/a_pb2,3)*(1.0-frac)*frac_s*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb2/a_pb3,3)*frac*(1.0-frac_s)*peb_map[i].mass_out[j];
	ring_mass_after+=pow(a_pb2/a_pb3,3)*(1.0-frac)*(1.0-frac_s)*peb_map[i].mass_out[j];
*/
	if(i<ring_num-1){
//	peb_map[i].mass_in[j]+=peb_map[i+1].fluxL[j]*dt1/dt0;
	}
	
	if(0 && sub_time1+2*dt1>dt0 && j>=9 && j<=11&& i<10 && i > 6) printf("PEB_DENS %d %d = %g ratio=%f\t%f\t%f dt=%f\n",i,j,peb_map[i].mass_in[j]/AREA,a_pb22/a_pb2,frac,frac_s,dt1);

}
	//for(j=0;j<1;j++){
        old_sigma=dust_budget[i].surf_dens;
     /*   dust_budget[i].mass_out-=(ring_mass_after-ring_mass_before);
        dust_budget[i].surf_dens-=(ring_mass_after-ring_mass_before)/AREA;
      */
	dust_budget[i].mass_out-=ring_mass_gain;
        dust_budget[i].surf_dens-=ring_mass_gain/AREA;
	if(0 && i<=i_lim1+2 && i >=i_lim1-2) printf("%d AREA=%e MASS_GAIN=%e RING=%e %e %e DUST=%e %e\n",i,AREA,ring_mass_gain,ring_mass_before,ring_sigma,ring_mass_before/AREA,dust_budget[i].surf_dens,dust_budget[i].mass_out);
	dust_budget[i].rho=dust_budget[i].rho*dust_budget[i].surf_dens/old_sigma;
        ratio_sigma=dust_budget[i].surf_dens/old_sigma;
        if(dust_budget[0].surf_dens<0.0) dust_budget[0].surf_dens=1e-10;

	if(dust_budget[i].surf_dens<0.0){
                printf("dust_surf_dens=%e\t%d\n",dust_budget[i].surf_dens,i);
                return -1.0;
        }
	//}
//        printf("HERE!!!\t%d\n",j);


        


        for(j=0;j<peb_size_num;j++){
		if(i<ring_num-1) peb_map[i].mass_in[j]+=peb_map[i+1].fluxL[j]*dt1/dt0;
		//if(i==0) peb_map[i].mass_in[j]+=peb_map[i].fluxL[j]*dt1/dt0; //wrong way
		if(i>0) peb_map[i].mass_in[j]+=peb_map[i-1].fluxR[j]*dt1/dt0;
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                //peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
								if (peb_map[i].surf_dens[j]<0.0) peb_map[i].surf_dens[j]=1e-30;
								peb_map[i].mass_out[j]=peb_map[i].surf_dens[j]*AREA;

        }
sub_time1+=dt1;
//	printf("SUB_TIME_AFT=%f\n",sub_time1);

}

}

return 1.0;
}


void dust_evolve(double dt0){

int i,i_new;
double vr_g, AREA,old_sigma,frac;
double factor=1.0;
dust_budget[ring_num-1].mass_out+=1.0*factor*mdot*MSUN*dt0*dust_gas;
for(i=ring_num-1;i>-1;i--){
	vr_g=vr_gas(dust_budget[i].rad)*1;	
	frac=vr_g*dt0*TUNIT/LUNIT/dust_budget[i].dr;
	i_new=i-1;
//	frac=0.0;
	if(i_new>-1) {
	dust_budget[i_new].mass_in+=frac*dust_budget[i].mass_out;
	dust_budget[i].mass_out-=frac*dust_budget[i].mass_out;
	}
	else dust_budget[i].mass_out-=frac*dust_budget[i].mass_out;

}

for(i=ring_num-1;i>-1;i--){
//	AREA=M_PI*((dust_budget[i].rad+size_ring/2.0)*(dust_budget[i].rad+size_ring/2.0)-(dust_budget[i].rad-size_ring/2.0)*(dust_budget[i].rad-size_ring/2.0))*LUNIT*LUNIT;	
	AREA=dust_budget[i].AREA;
	dust_budget[i].mass_out+=dust_budget[i].mass_in;
	dust_budget[i].mass_in=0.0;
	//if(i==0) dust_budget[i].mass_out=AREA*Sigma(dust_budget[i].rad)*dust_gas;
//	for(j=0;j<1;j++){
	old_sigma=dust_budget[i].surf_dens;
	dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
	dust_budget[i].rho=dust_budget[i].rho*dust_budget[i].surf_dens/old_sigma;
//}

}
}


void stokes_size(){
FILE* fp;
int i;
fp=fopen("max_size.txt","w");
for(i=0;i<ring_num;i++){
	fprintf(fp,"%g\t%g\n",peb_map[i].rad_med,2.25*mean_path(peb_map[i].rad_med));	
}
fclose(fp);
}

void tau_unity(){
FILE *fp;
int i;
double a1=0.1,a2=20000.0,a,vr0,tau=0.0;
fp=fopen("tau_unity.txt","w");
for(i=0;i<ring_num;i++){
	a1=0.1;
	a2=20000.0;
	tau=0.0;
	while(fabs(tau-1.0)>0.001){
	a=0.5*(a1+a2);
        vr0=vr_estimate(peb_map[i].rad_med,a,pp_vr_tau);
        tau=pp_vr_tau[1];
	if(tau>1.0) a2=a;
	else a1=a;
	}
	fprintf(fp,"%g\t%g\n",peb_map[i].rad_med,a);	
}
fclose(fp);
}


