//
//  drift.c
//  pebble_size_drift
//
//  Created by Xiao Hu on 5/16/15.
//
// velocity unit km/s

/*typedef struct PEBBLE{
	        double rad[10000];
		        double size[10000];
			        double time[10000];
				        double vr[10000];
} PEBBLE;*/

#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#define output_size "peb_size.txt"
//#define output_time "001alpha1cm100AU001AU1sun01acc.txt"
#define output_time "drift_test0.txt"
#define EFF 0.0
#define a_peb0 1.0 //radius in cm
//#define M_PI =3.141592654;
#define r_out 100.0
#define r_in 0.1
#define step0 0.001
#define acc_rate 1.0
#define int_m_star 1

#if int_m_star > 1
#define kk 0.57
#else
#define kk 0.8
#endif


double vr_tau[2]={0.0};
double tau_temp;


    
double a_peb(double r){
	return a_peb0*1.0;//*np.exp((r_out-r)/5);
}


double tau_fric(double r, double a_pb){
	    alpha = alpha_func(r);
    double t_eps;
    t_eps=0.0178*a_pb*rho_peb*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,0.8)*\
    pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,-0.6)*pow(r,0.6);
    
    // if (a_pb < 2.25*mean_path(r)) {
    return  0.0178*a_pb*rho_peb*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,0.8)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,(-0.6))*pow(r,0.6);
    //}
    //else{
    //  return  (0.2222*a_pb*sound_sp(r)*100000.0/viscosity(r))*t_eps;
    //}
}

double tau_test(double r, double a_pb){
    //return   w_K(r)*rho_peb*a_pb/(density(r)*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*100000.0);
    return     w_K(r)*2.0*rho_peb*a_pb*a_pb/(9.0*viscosity(r)*density(r));
    
    
}

double tau_fric0(double r, double a_pb){
    double t_eps;
		alpha = alpha_func(r);
    t_eps=0.0178*a_pb*rho_peb*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,0.8)*\
    pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,-0.6)*pow(r,0.6);
    
    if ( a_pb < 2.25*mean_path(r)) {
        //return  0.0178*a_pb*rho_peb*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,0.8)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,(-0.6))*pow(r,0.6);
        return w_K(r)*rho_peb*a_pb/(density(r)*sqrt(8.0/gamma0/M_PI)*sound_sp(r));
        
    }
    else{
        //return  (0.2222*a_pb*sqrt(8.0/gamma0/M_PI)*sound_sp(r)/viscosity(r))*t_eps*100000.0;
        return 2.0*w_K(r)*rho_peb*a_pb*a_pb/(9.0*viscosity(r)*density(r));
    }
}

double v_r00(double r,double a_pb){
	    alpha = alpha_func(r);
    //printf("tau=%f\n",tau_fric(r));
    return 100000.0*0.108*(1.0/(tau_fric0(r,a_pb)+1.0/tau_fric0(r,a_pb)))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
    
    
}

double v_r0(double r, double a_pb){
    // printf("tau=%f\n",tau_fric(r));
    //return 0.108*(1.0/(tau_fric(r)+1.0/tau_fric(r)))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
    double tau;
    tau=w_K(r)*rho_peb*a_pb/(density(r)*sqrt(8.0/gamma0/M_PI)*sound_sp(r));
    //return 0.108*(1.0/(tau+1.0/tau))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
    return yeta(r)/(tau+1.0/tau)*v_K(r);
    
    
    
    
}

double v_r1(double r, double a_peb){
    double tau;
    //tau=(0.2222*a_pb*sqrt(8.0/gamma0/M_PI)*sound_sp(r)/viscosity(r))*tau_fric(r)*100000.0;
    tau=w_K(r)*2.0*rho_peb*a_peb*a_peb/(9.0*viscosity(r)*density(r));
    //tau=tau_fric0(r);
    //printf("tau=%f\n",tau);
    //printf("tau=%1.20f\n",tau);
    // return 0.108*(1.0/(tau+1.0/tau))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
    return yeta(r)/(tau+1.0/tau)*v_K(r);
    
}

double v_peb_gas(double r, double a_pb){//in cm/s
    return 0.5*tau_fric(r,a_pb)*v_r0(r,a_pb)*v_K(r)*pow(10.0,10.0)/(w_K(r)*r*AU_km*pow(10,5));
}

double Reynolds(double r, double vr, double tau, double a_pb){
    double v_theta;
    v_theta=0.5*tau*vr;
    //return 2*a_pb*sqrt(v_theta*v_theta+vr*vr*1e10)/viscosity(r);
    return 2*a_pb*sqrt(v_theta*v_theta+vr*vr)/viscosity(r);
    
}

double *v_r2(double r, double a_pb){
    double v1,v2,f0,f1,tau0,tau,tau1,t1,t2,re0;
    int j,i;
    v1=0.07*100000;
    v2=0.044556*100000;
    t1=0.3;
    if (Reynolds(r,v_r1(r,a_pb),tau_fric0(r,a_pb),a_pb)-1.0<0.1) t2=tau_fric0(r,a_pb);
    else t2=tau_temp;
    j=0;
    i=0;
    while(fabs(t2-t1)>fabs(t1)*0.000000000000001){
  //      if (j==0) t1=t2;
    //    else t1=0.5*t1+0.5*t2;
	    t1=t2;
        //t2=3.095783;
        // t1=0.038;
        i=0;
        while(fabs(v2-v1)>fabs(v1)*0.0000001){// && fabs(t2-t1)>fabs(t1)*0.0001){
            //if(i==0) v1=v2;
            //else v1=0.5*v1+0.5*v2;
		v1=v2;
            //t1=t2;
            // t1=0;
            re0=2.0*a_pb*sqrt(1.0+0.25*t1*t1)*v1/viscosity(r);
            //printf("re0=%f\n",re0);
            tau0=w_K(r)*rho_peb*pow(a_pb,1.6)/(9.0*density(r)*pow(viscosity(r),0.6))*pow(sqrt(1.0+0.25*t1*t1),-0.4)*pow(2.0,0.6);
            //tau0=w_K(r)*rho_peb*pow(a_pb,1.6)/(72.0*density(r))*pow(sqrt(1.0+0.0*t1*t1)/viscosity(r),0.6)*1.51571656651;
            tau=tau0*pow(v1,-0.4);
            tau1=-0.4*tau0*pow(v1,-1.4);//derivative of tau over v1
            f0=yeta(r)/(tau+1.0/tau)-v1/v_K(r);
            f1=-1.0*yeta(r)/(tau+1.0/tau)/(tau+1.0/tau)*(tau1-tau1/tau/tau)-1.0/v_K(r);
            //t2=tau;
            v2=v1-f0/f1;
            if (v2<0.0001 && v2>0.0) {
                v2+=0.0001;
            }
            else if(v2<0.0){
                v2*=-1.0;        }
            //v1=0.035277;
            //v2=0.035;
            i++;
            //printf("tau0=%1.20f\n",tau);
//	    printf("V_R2 i=%d\tj=%d\n",i,j);
		t1=tau;
            
        }
        t2=tau;
        j++;
    }
    vr_tau[0]=v2;
    vr_tau[1]=tau0*pow(v2,-0.4);
//    printf("V_R2 i=%d\tj=%d\n",i,j);
    return vr_tau;
}

double *v_r3(double r, double a_pb){
    //coeff=w_K(r)*rho_peb*a_pb/(0.44*density(r));
    // printf("%f\n",(coeff*(1.0*yeta(r)*v_K(r)+coeff)));
    //printf("tau_temp=%f\n",tau_temp);
    //  printf("gua3");
    double v1,v2,tau0,t1,t2,ft,ft1,fx,fx1,B;
    int i,j=0;
    v1=0.12*10000000;;
    v2=0.125*10000000;
    t1=50.0;
    //if (Reynolds(r,v_r2(r)[0],v_r2(r)[1])-800.0<100 || 1) t2=v_r2(r)[1];
    t2=tau_temp;
    //t2=2.06;
    //printf("vr3aaaaaa %f \n",t2);
    //printf("%f\t%e\t%e\t%e\t%e\n",r,w_K(r),v_K(r),density(r),yeta(r));
    while(fabs(t2-t1)>fabs(t1)*0.00001){
        //if (j==0) t1=t2;
        //else if(fabs(t2-t1)>0.05*fabs(t1)) t1=0.01*t2+0.99*t1;
        //else t1=t2*1.0+t1*0.0;
        t1=0.0*t1+1.0*t2;
        //printf("t1=%f\n",t1);
        //while(fabs(v2-v1)>fabs(v1)*0.0001){// && fabs(t2-t1)>fabs(t1)*0.0001){
        // v1=v2;
        //t1=t2;
        // t1=0;
        //t1=247;
        // printf("t1=%f\n",t1);
        ft=w_K(r)*4.0/3.0*rho_peb*a_pb/(0.5*24*pow(800.0,-0.6)*density(r)*sqrt(1.0+0.25*t1*t1));
        ft1=-0.25*t1*ft/(1+0.25*t1*t1);
        B=yeta(r)*v_K(r);
        fx=ft/sqrt(fabs(-1.0*ft*ft+ft*B))-t1;
        if (-1.0*ft*ft+ft*yeta(r)*v_K(r) >=0){
            fx1=-0.5*pow((B*ft-ft*ft),-1.5)*(-2.0*ft*ft1+B*ft1)*ft+ft1/sqrt(B*ft-ft*ft)-1.0;
        }
        else fx1=-0.5*pow((-1.0*B*ft+ft*ft),-1.5)*(2.0*ft*ft1-B*ft1)*ft+ft1/sqrt(-1.0*B*ft+ft*ft)-1.0;
        t2=t1-fx/fx1;
        tau0=ft;
        //tau0=w_K(r)*4.0/3.0*rho_peb*a_pb/(0.5*24*pow(800.0,-0.6)*density(r)*sqrt(1.0+0.25*t1*t1));
        //tau0=w_K(r)*rho_peb*pow(a_pb,1.6)/(72.0*density(r))*pow(sqrt(1.0+0.0*t1*t1)/viscosity(r),0.6)*1.51571656651;
        
        //tau=tau0/v1;
        //tau1=-1.0*tau0/v1/v1;
        //f0=yeta(r)/(tau+1.0/tau)-v1/v_K(r)/100000;
        //f1=-1.0*yeta(r)/(tau+1.0/tau)/(tau+1.0/tau)*(tau1-tau1/tau/tau)-1.0/v_K(r)/100000;
        //t2=tau;
        //v2=v1-f0/f1;
        //if (v2<0.000000000000000000001) {
        //     v2+=0.000000000000000001;
        // }
        // printf("tau=%1.20f\n",v1);
        
        //}
        j++;
        //v2=sqrt(fabs(yeta(r)*tau0*v_K(r)*100000.0-tau0*tau0));
        v2=sqrt(fabs(-1.0*ft*ft+ft*B));
        v2=yeta(r)*v_K(r)/(t1+1.0/t1);
        //v2=0.123968*100000.0;
        //t2=tau0/v2;
        // printf("r=%f\tcoeff=%f\t tau=%f\t %f\n",r,t2,fx,fx1);
        
    }
    fx=10.0;
    if (fabs(t2-tau_temp)>0.1*tau_temp) {
        double t_temp;//[1000]={0};
        //printf("we tried \t%f\n",tau_temp);
        i=0;
        while (fabs(fx)>1.0 && i<10000){
            t_temp=tau_temp*(1.0+0.4*(i-5000)/10000);
            t1=t_temp;
            ft=w_K(r)*4.0/3.0*rho_peb*a_pb/(0.5*24*pow(800.0,-0.6)*density(r)*sqrt(1.0+0.25*t1*t1));
            ft1=-0.25*t1*ft/(1+0.25*t1*t1);
            B=yeta(r)*v_K(r);
            fx=ft/sqrt(fabs(-1.0*ft*ft+ft*B))-t1;
            i++;
        }
        
    }
    t2=t1;
    v2=yeta(r)*v_K(r)/(t1+1.0/t1);
    //printf("r=%f\t tau0=%f\t tau=%f\t fx=%f\n",r,tau0,t2,fx);
    vr_tau[0]=v2;//sqrt(coeff*(yeta(r)*v_K(r)+coeff));
    vr_tau[1]=t2;//w_K(r)*rho_peb*a_pb/1.32/density(r)/vr_tau[0];
    // printf("vr3bbb \n");
    return vr_tau;
    
    
}

/*double v_r(double r){
 double v1,v2,f0,f1,tau0,tau,tau1,coeff;
 //printf("%f\n",mean_path(r));
 v1=200;
 v2=100;
 tau0=w_K(r)*rho_peb*3.0*pow(a_pb,1.6)/(72.0*density(r))*pow(viscosity(r),-0.6)*1.51571656651;
 if (a_pb < 2.25*mean_path(r)|| Reynolds(r)<=1.0 ) {
 //printf("gua1");
 return 0.108*(1.0/(tau_fric(r)+1.0/tau_fric(r)))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
 
 }
 
 else if (Reynolds(r) < 800.0 && Reynolds(r)>1.0) {
 while(fabs(v2-v1)>fabs(v1)*0.0001){
 v1=v2;
 tau=tau0*pow(v1,-0.4);
 tau1=-0.4*tau0*pow(v1,-1.4);
 f0=yeta(r)/(tau0*pow(v1,-0.4)+1.0/(tau0*pow(v1,-0.4)))-v1/v_K(r);
 f1=-1.0*yeta(r)/(tau+1.0/tau)/(tau+1.0/tau)*(tau1-tau1/tau/tau)-1.0/v_K(r);
 v2=v1-f0/f1;
 }
 //printf("v2=%f\n",v2);
 //printf("%f\n",yeta(r));
 //printf("%f\n",tau0);
 //printf("%f\n",v_K(r));
 //printf("%f\n",r);
 return v2;
 }
 else{
 coeff=w_K(r)*rho_peb*a_pb/(0.44*density(r));
 printf("%f\n",(coeff*(1.0*yeta(r)*v_K(r)+coeff)));
 printf("%f\n",yeta(r));
 printf("gua3");
 return sqrt(coeff*(yeta(r)*v_K(r)+coeff));
 
 }
 
 }*/
//    return r**-2


double v_r1_test(double r, double a_pb){
    double tau;
    //tau=(0.2222*a_pb*sqrt(8.0/gamma0/M_PI)*sound_sp(r)/viscosity(r))*tau_fric(r)*100000.0;
    tau=w_K(r)*2.0*rho_peb*a_pb*a_pb/(9.0*viscosity(r)*density(r));
    //tau=tau_fric0(r);
    //printf("tau=%f\n",tau);
    //printf("tau=%1.20f\n",tau);
    // return 0.108*(1.0/(tau+1.0/tau))*pow(gamma0,0.8)*pow(opa,0.2)*pow(alpha,-0.2)*pow(m_star,-0.2)*pow((1-sqrt(r_star/r))*acc_rate,0.4)*pow(r,-0.4);
    return yeta(r)/(tau+1.0/tau)*v_K(r);

}

double vr_estimate( double r, double a_pb, double *p_vr_tau)
{
	double drag_coeff,drag_force,Re,Re2,tau_fric,vr,vt,output;
	int judge=0;
	Re2=2.0*a_pb*(v_K(r)-vt_gas(r))/viscosity(r);
	do{
	Re=Re2;
	if(a_pb<2.25*mean_path(r)){
		drag_coeff=8.0/3.0*sqrt(8.0/gamma0/M_PI)*sound_sp(r)/(v_K(r)-vt_gas(r));
		judge=1;
	}
	else if(Re<=1.0){
		drag_coeff=24.0/Re;
		judge=2;
	}
	else if(Re>1.0 && Re < 800.0){
		drag_coeff=24.0*pow(Re,-0.6);
		judge=3;
	}
	else{
		drag_coeff=0.44;
		judge=4;
	}
	drag_force=0.5*drag_coeff*M_PI*a_pb*a_pb*density(r)*(v_K(r)-vt_gas(r))*(v_K(r)-vt_gas(r));
//	return drag_force*r*LUNIT/(rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*v_K(r));
	tau_fric=rho_peb*4.0/3.0*M_PI*pow(a_pb,3)*(v_K(r)-vt_gas(r))/drag_force*w_K(r);
//	tau=w_K(r)*2.0*rho_peb*a_pb*a_pb/(9.0*viscosity(r)*density(r));
//	vr=yeta(r)/(tau+1.0/tau)*v_K(r);
	vr=yeta(r)*v_K(r)/(tau_fric+1.0/tau_fric);
	vt=0.5*tau_fric*vr;
	Re2=2.0*a_pb*sqrt(vr*vr+vt*vt)/viscosity(r);
	}while(fabs(Re-Re2)>0.0001*fabs(Re2) && 0);
	//vr=v_r1_test(r,a_pb);
	output=yeta(r)*v_K(r)/(tau_fric+1.0/tau_fric);
	output=vr;
	//if((a_peb>2.25*mean_path(r) && Re<=1.0) && fabs(tau-tau_fric)>0.01*tau) printf("DIFF tau_fric=%f tau=%f\n",tau_fric,tau);
        //if(judge==2 || 1) printf("DIFF tau_fric=%f tau=%f\t%d\t%f\n",output,vr,judge,Re);

	p_vr_tau[0]=output;
	p_vr_tau[1]=tau_fric;
//	printf("r=%g\t tau_fric=%g\t Omega_K=%g\t v_K=%g\n",r,tau_fric,v_K(r)/r/LUNIT,v_K(1.0));
	return output;

}
	
void vr_iter(double r, double a_peb, double vr1, double vt1, double *vr2, double *vt2, double *tau_fric){
		double drag_coeff,drag_force,Re;
                Re=2.0*a_peb*sqrt(vr1*vr1+vt1*vt1)/viscosity(r);
		//power = Re<=1.0?-1.0:-0.6;
		if(Re <=1.0) drag_coeff=24.0/Re;
		else if( Re> 1.0 && Re <=800) drag_coeff=24.0*pow(Re,-0.6);
		else drag_coeff=0.44;
                //printf("%f %f\t ",Re,power);
                //drag_coeff=24.0*pow(Re,power);
                drag_force=0.5*drag_coeff*M_PI*a_peb*a_peb*density(r)*(vr1*vr1+vt1*vt1);
                *tau_fric=rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*sqrt(vr1*vr1+vt1*vt1)/drag_force*w_K(r);
                *vr2=yeta(r)*v_K(r)/(*tau_fric+1.0/(*tau_fric));
                *vt2=0.5*(*tau_fric)*(*vr2);
  }

double drift_vr(double r,double a_peb, double *p_vr_tau)
{
	double vr,vt,vr1,vr2,vt1,vt2;
        double drag_force,Re,tau_fric;
	int count;
	vr=vr_estimate(r,a_peb,p_vr_tau);
	tau_fric=p_vr_tau[1];
	vt=0.5*tau_fric*vr;
	vr2=vr;
	vt2=vt;
	vr1=vr;
	vt1=vt;
        Re=2.0*a_peb*sqrt(vr*vr+vt*vt)/viscosity(r);
	if(0 || a_peb<2.25*mean_path(r) ){
        p_vr_tau[0]=vr;
        p_vr_tau[1]=tau_fric;
	return vr;
	}
	else if(Re<=1.0 && 0){
	count=0;
        do{
                vr1=0.1*vr2+0.9*vr1;
                vt1=0.1*vt2+0.9*vt1;
                vr_iter(r,a_peb,vr1,vt1,&vr2,&vt2,&tau_fric);
                count++;
        }while(fabs(vr1-vr2)>fabs(1e-4*vr1));
	p_vr_tau[0]=vr2;
        p_vr_tau[1]=rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*sqrt(vr2*vr2+vt2*vt2)/drag_force*w_K(r);
        return vr2;
	}
	
	else if(1 || Re<800.0){
	count=0;
        do{
		vr1=0.5*vr2+0.5*vr1;
                vt1=0.5*vt2+0.5*vt1;
		vr_iter(r,a_peb,vr1,vt1,&vr2,&vt2,&tau_fric);
		count++;	
	//	if(count>1000)	printf("vr_iter=%d\t %f\t %f\n",count,vr1,vr2);
		if(count>10000) break;

	}while(fabs(vr1-vr2)>1e-4*vr2);
	//printf("vr_iter=%d\n",count);
	p_vr_tau[0]=vr2;
        p_vr_tau[1]=tau_fric;
	return vr2;
	}
	else {
	printf("ERROR!\n");
	return -1.0;
	}
		

}


double drift_vr_test(double r,double a_peb, double *p_vr_tau)
{
        double vr,vt,vr1,vr2,vt1,vt2,v;
        double drag_coeff,drag_force,Re,tau_fric;
        int count,i;
	FILE *fp;
        vr=vr_estimate(r,a_peb,p_vr_tau);
        tau_fric=p_vr_tau[1];
        vt=0.5*tau_fric*vr;
        vr2=vr;
        vt2=vt;
        vr1=vr;
        vt1=vt;
        Re=2.0*a_peb*sqrt(vr*vr+vt*vt)/viscosity(r);
	fp=fopen("vr_test.txt","w");
	for(i=1;i<40000;i++){
		v=i*1.0;
		Re=2.0*a_peb*v/viscosity(r);
		drag_coeff=24.0*pow(Re,-0.6);
                drag_force=0.5*drag_coeff*M_PI*a_peb*a_peb*density(r)*v*v;
                tau_fric=rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*v/drag_force*w_K(r);
                vr=yeta(r)*v_K(r)/(tau_fric+1.0/tau_fric);
                vt=0.5*tau_fric*vr;
		fprintf(fp,"%e\t%e\n",v/sqrt(1.0+0.25*tau_fric*tau_fric),v/sqrt(1.0+0.25*tau_fric*tau_fric)-vr);	
	}
	fclose(fp);
	
        if(0 || a_peb<2.25*mean_path(r) || Re<=1.0){
        p_vr_tau[0]=vr;
        p_vr_tau[1]=tau_fric;
        return vr;
        }
        else if(Re<800.0){
        count=0;
        do{
                vr1=0.1*vr2+0.9*vr1;
                vt1=0.1*vt2+0.9*vt1;
                Re=2.0*a_peb*sqrt(vr1*vr1+vt1*vt1)/viscosity(r);
                printf("%e ",vr1);
                drag_coeff=24.0*pow(Re,-1.0);
                drag_force=0.5*drag_coeff*M_PI*a_peb*a_peb*density(r)*(vr1*vr1+vt1*vt1);
                tau_fric=rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*sqrt(vr1*vr1+vt1*vt1)/drag_force*w_K(r);
                vr2=yeta(r)*v_K(r)/(tau_fric+1.0/tau_fric);
                vt2=0.5*tau_fric*vr2;
                count++;
        }while(fabs(vr1-vr2)>fabs(0.001*vr1));
        printf("%d\n",count);
        p_vr_tau[0]=vr2;
        p_vr_tau[1]=rho_peb*4.0/3.0*M_PI*pow(a_peb,3)*sqrt(vr2*vr2+vt2*vt2)/drag_force*w_K(r);
        return vr2;
        }
	
    return 0.0;


}






int drift(double r_start, double a_pebble, double coag_eff)
{
    
    int i=0,j=0,k=0,l=0;
    FILE *fp_vr,*fp_drt,*fp_size;
    fp_vr=fopen("drift_velocity", "w");
    fp_drt=fopen(output_time,"w");
    fp_size=fopen(output_size,"w");
    
    double x0,x1,x,x_cut,x_stop,y,Re1,Re2,vr0,vr1,vr2,a_pb1,a_pb2,tau,vol_plus;
    double k1,k2,k3,k4,step,sum1=0.0;
    x0=r_start;
    
    x1=r_in;x=x0;y=0.0;
    //coag_eff=EFF;
    step=-1.0*step0;
    a_pb2=a_pebble;
    while (x>x1) {
        if(tau_temp <1e8 || 1)
        {
            j=0;
            a_pb1=a_pb2;
            //printf("pass0 \n");
            Re1=Reynolds(x,v_r0(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
            Re2=0.0;
            //printf("Re1=%f\n",Re1);
            if (Re1 >2.0) {
                Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
                // Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
                //Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
                //tmp1=v_r2(x,a_pb1)[0];
                //  tmp2=v_r2(x,a_pb1)[1];
                //   Re2=Reynolds(x, tmp1,tmp2,a_pb1);
                //Re1=v_r2(x,a_pb1)[0];
                // printf("Re2=%f\n",Re2);
                // printf("Stokes2b Re=%f\tVr=%f\t TAU=%f\t Yeta=%f\t x=%f\t apb=%f\n",Re2,v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],yeta(x),x,a_pb1);
                //tmp1=v_r2(x,a_pb1)[0];tmp2=v_r2(x,a_pb1)[1];
                
                //  Re2=Reynolds(x, tmp1,tmp2,a_pb1);
                //  printf("Stokes2a Re=%f\tVr=%f\t TAU=%f\t Yeta=%f\t x=%f\t apb=%f\n",Re2,v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],yeta(x),x,a_pb1);
                
                
            }
            while (fabs(Re1-Re2)>0.001*fabs(Re1)){
                Re1=Re2;
                if (a_pb1 < 2.25*mean_path(x)) {
                    Re2=Reynolds(x, v_r0(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
                    //printf("Epstein Re=%f\tVr=%f TAU=%f\t Yeta=%f\n",Re2,v_r0(x,a_pb1),tau_fric0(x,a_pb1),yeta(x));
                }
                else if (Re1 <= 1.0 && 1){
                    Re2=Reynolds(x, v_r1(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
                    //printf("Stokes1 Re=%f\tVr=%f\t TAU=%f\t Yeta=%f\n",Re2,v_r1(x,a_pb1),tau_fric0(x,a_pb1),yeta(x));
                    
                }
                else if (Re1 <=800.0 && Re1>1.0 && 1 ) {
                    Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
                    
                    //printf("Stokes2 Re=%f\tVr=%f\t TAU=%f\t Yeta=%f\t x=%f\t apb=%f\n",Re2,v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],yeta(x),x,a_pb1);
                }
                else if (Re1 > 800.0 && 1 ){
                    Re2=Reynolds(x, v_r3(x,a_pb1)[0],v_r3(x,a_pb1)[1],a_pb1);
                    //printf("Stokes3 Re=%f\tVr=%g\t TAU=%f\n",Re2,v_r3(x,a_pb1)[0],v_r3(x,a_pb1)[1]);
                }
                j++;
                if (j>11) {
                    //break;
                    // printf("Re=%f \trad=%f\n",Re2,x);
                }
                
            }
            //printf("Re=%f\n",Re2);
            
            if (a_pb1 < 2.25*mean_path(x)) {
                vr0=v_r0(x,a_pb1);
                vr1=v_r0(x+0.5*step,a_pb1);
                vr2=v_r0(x+step,a_pb1);
            }
            else if (Re2<=1.0 && 1){
                vr0=v_r1(x,a_pb1);
                vr1=v_r1(x+0.5*step,a_pb1);
                vr2=v_r1(x+step,a_pb1);
            }
            else if (Re2 < 800.0 && Re2>1.0 && 1 ) {
                vr0=v_r2(x,a_pb1)[0];
                vr1=v_r2(x+0.5*step,a_pb1)[0];
                vr2=v_r2(x+step,a_pb1)[0];
            }
            else if (Re2 >800.0 && 1 ){
                vr0=v_r3(x,a_pb1)[0];
                vr1=v_r3(x+0.5*step,a_pb1)[0];
                vr2=v_r3(x+step,a_pb1)[0];
            }
            if(fabs(vr0-vr2)>0.01*fabs(vr0)) {
                //if (0) {
                step=step*0.8;
                //step=step*0.5;
                if (fabs(step) < 0.0001) step=-0.0001;
                //printf("change step to %1.12f\t vr accuracy=%f\n",step,fabs(vr0-vr2)/fabs(vr0));
                //printf("change step to %1.12f\n",step);
            }
            if (a_pb1 < 2.25*mean_path(x) && 1 ) {
                vr0=v_r0(x,a_pb1);
                vr1=v_r0(x+0.5*step,a_pb1);
                vr2=v_r0(x+step,a_pb1);
                tau=tau_fric0(x,a_pb1);
            }
            else if (Re2<=1.0 && 1){
                vr0=v_r1(x,a_pb1);
                vr1=v_r1(x+0.5*step,a_pb1);
                vr2=v_r1(x+step,a_pb1);
                tau=tau_fric0(x,a_pb1);
            }
            else if (Re2 < 800.0 && Re2>1.0 || 0 ) {
                vr0=v_r2(x,a_pb1)[0];
                vr1=v_r2(x+0.5*step,a_pb1)[0];
                vr2=v_r2(x+step,a_pb1)[0];
                tau=v_r2(x,a_pb1)[1];
            }
            else if (Re2 >800.0 && 1 ){
                vr0=v_r3(x,a_pb1)[0];
                vr1=v_r3(x+0.5*step,a_pb1)[0];
                vr2=v_r3(x+step,a_pb1)[0];
                tau=v_r3(x,a_pb1)[1];
            }
            tau_temp=tau;
            k1=1.0/vr0;
            k2=1.0/vr1;
            k3=1.0/vr1;
            k4=1.0/vr2;
            y=y+step*(k1+2*k2+2*k3+k4)/6.0;
            fprintf(fp_vr, "%f\t%f\n",x,vr0);
            //fprintf(fp_drt, "%f\t%f\t%f\t%f\n",x,-1.0*sum1*AU_km/TUNIT,a_pb1,coag_eff);
            fprintf(fp_drt, "%f\t%f\n",x,-1.0*sum1*LUNIT/TUNIT);
            fprintf(fp_size,"%f\t%f\n",x,a_pb1);
            if (1) {
                //printf("r=%0.20f\t tau=%f\t f_tau=%f\t t_drift=%f\t vr=%f\tpeb_size=%f\t Step=%0.20f\n",x,tau,1.0/(tau+1.0/tau),43.9*(tau+1.0/tau)*pow(acc_rate,-0.4)*pow(alpha,0.2)*pow(x,1.4), vr0,0.108/(tau+1.0/tau)*pow(acc_rate,0.4)*pow(alpha,-0.2)*pow(x,-0.4),step);
            }
            sum1=sum1+step*k2;
            x+=step;
            if (a_pb1>9.0/4.0*mean_path(x)) {
                l++;
                if(l==1) x_stop=x;
                coag_eff=coag_eff*exp((x-x_stop)/0.001);
                coag_eff=0.000;
                
            }
            vol_plus=-1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*step*(k1+2*k2+2*k3+k4)/6.0*LUNIT;
            if (0) vol_plus=0.0;
            a_pb2=pow(((vol_plus*coag_eff*density(x)/rho_peb+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
            /*if (i%100==0) {
             printf("%0.20f\t peb_size=%f\t Step=%0.20f\n",x,a_pb2,step);
             }*/
            //printf("end \n");
        }
        else{
            k++;
            if(k==1) x_cut=x;
            y=y+step/vr0;
            x+=step;
            fprintf(fp_vr, "%f\t%f\n",x_cut,vr0);
            //fprintf(fp_drt, "%f\t%f\t%f\t%f\n",x_cut,-1.0*sum1*AU_km/TUNIT,a_pb1,coag_eff);
            fprintf(fp_drt, "%f\t%f\t%f\t%f\n",x_cut,-1.0*sum1*AU_km/TUNIT,a_pb1,coag_eff);
            fprintf(fp_size,"%f\t%f\n",x_cut,a_pb1);
            
        }
        
    }
    fclose(fp_vr);
    fclose(fp_drt);
    fclose(fp_size);
    //printf("%transition=%f\n",x_stop);
    
    
    sum1=sum1*LUNIT/TUNIT;
    
    //printf("%0.10f\n",sum1);
    //printf("Hello, World!\n");
    
    x0=r_out;x1=r_in;x=x0;y=0.0;
    step=-1.0*step0;
    sum1=0.0;
    while (x>x1) {
        if(fabs(v_r00(x,a_peb(x))-v_r00(x+step,a_peb(x)))>0.001*fabs(v_r00(x,a_peb(x)))) {
            //if (0) {
            step=step*0.5;
            if (fabs(step) < 0.000001) step=-0.000001;
            //printf("change step to %1.12f\t vr accuracy=%f\n",step,fabs(v_r00(x,a_peb(x))-v_r00(x+step,a_peb(x)))/fabs(v_r00(x,a_peb(x))));
        }
        k1=1.0/v_r00(x,a_peb(x));
        k2=1.0/v_r00(x+0.5*step,a_peb(x));
        k3=1.0/v_r00(x+0.5*step,a_peb(x));
        k4=1.0/v_r00(x+step,a_peb(x));
        y=y+step*(k1+2*k2+2*k3+k4)/6.0;
        sum1=sum1+step*k2;
        x+=step;
        if (i%100==0) {
            //printf("%f\r",(x0-x)/(x0-x1)*100.0);
        }
        
    }
    sum1=sum1*AU_km/TUNIT;
    printf("%0.10f\n", sum1);
    printf("%0.10f\n", tau_fric0(1.0,a_pb1));
    printf("%0.10f\n", v_K(0.01516));
    printf("%0.10f\n", yeta(0.01516));
    
    
    
    return 0;
}



/*int drift_t(PEBBLE *pp, double coag_eff, int num)
{

    int i=0,j=0,k=0,l=0,ll=0,zone_num1,zone_num2;
    FILE *fp_vr,*fp_peb,*fp_size;
    char output_peb[256];
    fp_vr=fopen("drift_velocity", "w");
    sprintf(output_peb,"pebble_num_%d.txt",num);
    fp_peb=fopen(output_peb,"w");
    fp_size=fopen(output_size,"w");
    double x0,x1,x,x_cut,x_stop,y,Re1,Re2,vr0,vr1,vr2,a_pb1,a_pb2,tau,vol_plus,tmp1,tmp2,time_tot=0.0,ring_area,ring_area0,ring_area_all,ring_mass,budget,massplus;
    double k1,k2,k3,k4,step,sum1=0.0;
    x0=pp->rad[0];
    
    x1=r_in;x=x0;y=0.0;
    //coag_eff=EFF;
    step=-1.0*step0;
    a_pb2=pp->size[0];
    while (x>x1) {
        if(tau_temp <1e8 || 1)
        {
            j=0;
            a_pb1=a_pb2;
            Re1=Reynolds(x,v_r0(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
            Re2=0.0;
            if (Re1 >2.0) {
                Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
            }
            while (fabs(Re1-Re2)>0.001*fabs(Re1)){
                Re1=Re2;
                if (a_pb1 < 2.25*mean_path(x)) {
                    Re2=Reynolds(x, v_r0(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
                }
                else if (Re1 <= 1.0 && 1){
                    Re2=Reynolds(x, v_r1(x,a_pb1),tau_fric0(x,a_pb1),a_pb1);
                    
                }
                else if (Re1 <=800.0 && Re1>1.0 && 1 ) {
                    Re2=Reynolds(x, v_r2(x,a_pb1)[0],v_r2(x,a_pb1)[1],a_pb1);
                    
                }
                else if (Re1 > 800.0 && 1 ){
                    Re2=Reynolds(x, v_r3(x,a_pb1)[0],v_r3(x,a_pb1)[1],a_pb1);
                }
                j++;
                if (j>11) {
                }
                
            }
            
            if (a_pb1 < 2.25*mean_path(x)) {
                vr0=v_r0(x,a_pb1);
                vr1=v_r0(x+0.5*step,a_pb1);
                vr2=v_r0(x+step,a_pb1);
            }
            else if (Re2<=1.0 && 1){
                vr0=v_r1(x,a_pb1);
                vr1=v_r1(x+0.5*step,a_pb1);
                vr2=v_r1(x+step,a_pb1);
            }
            else if (Re2 < 800.0 && Re2>1.0 && 1 ) {
                vr0=v_r2(x,a_pb1)[0];
                vr1=v_r2(x+0.5*step,a_pb1)[0];
                vr2=v_r2(x+step,a_pb1)[0];
            }
            else if (Re2 >800.0 && 1 ){
                vr0=v_r3(x,a_pb1)[0];
                vr1=v_r3(x+0.5*step,a_pb1)[0];
                vr2=v_r3(x+step,a_pb1)[0];
            }
            if(fabs(vr0-vr2)>0.01*fabs(vr0)) {
                step=step*0.8;
                if (fabs(step) < 0.0001) step=-0.0001;
            }
            if (a_pb1 < 2.25*mean_path(x) && 1 ) {
                vr0=v_r0(x,a_pb1);
                vr1=v_r0(x+0.5*step,a_pb1);
                vr2=v_r0(x+step,a_pb1);
                tau=tau_fric0(x,a_pb1);
            }
            else if (Re2<=1.0 && 1){
                vr0=v_r1(x,a_pb1);
                vr1=v_r1(x+0.5*step,a_pb1);
                vr2=v_r1(x+step,a_pb1);
                tau=tau_fric0(x,a_pb1);
            }
            else if (Re2 < 800.0 && Re2>1.0 || 0 ) {
                vr0=v_r2(x,a_pb1)[0];
                vr1=v_r2(x+0.5*step,a_pb1)[0];
                vr2=v_r2(x+step,a_pb1)[0];
                tau=v_r2(x,a_pb1)[1];
            }
            else if (Re2 >800.0 && 1 ){
                vr0=v_r3(x,a_pb1)[0];
                vr1=v_r3(x+0.5*step,a_pb1)[0];
                vr2=v_r3(x+step,a_pb1)[0];
                tau=v_r3(x,a_pb1)[1];
            }
            tau_temp=tau;
            k1=1.0/vr0;
            k2=1.0/vr1;
            k3=1.0/vr1;
            k4=1.0/vr2;
            y=y+step*(k1+2*k2+2*k3+k4)/6.0;
	    pp->vr[ll]=vr0;
	    pp->size[ll]=a_pb1;
	    pp->rad[ll]=x;
	    pp->time[ll]=time_tot;
	    pp->mass[ll]=pow(pp->size[ll]/pp->size[0],3)*pp->mass[0];
	    zone_num1=floor(x/dust_budget[0].dr);

	    x+=-vr0*dt/AU_km*TUNIT;

	    zone_num2=floor(x/dust_budget[0].dr);
	    if(zone_num2==zone_num1) {
	    	ring_area=(AU_km*100000.0)*(AU_km*100000.0)*M_PI*((x+vr0*dt/AU_km*TUNIT)*(x+vr0*dt/AU_km*TUNIT)-x*x);
	    	budget=dust_budget[zone_num1].surf_dens*ring_area;
	    }
	    else{
		    budget=0.0;
		    for(i=zone_num1-1;i>zone_num2;i--){
			    ring_area=(dust_budget[i].rad+0.125)*(dust_budget[i].rad+0.125)-(dust_budget[i].rad-0.125)*(dust_budget[i].rad-0.125);
			    ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
			    budget+=dust_budget[i].surf_dens*ring_area;
		    }
		    ring_area=(x+vr0*dt/AU_km*TUNIT)*(x+vr0*dt/AU_km*TUNIT)-(dust_budget[zone_num1].rad-0.125)*(dust_budget[zone_num1].rad-0.125);
		    ring_area=ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		    budget+=dust_budget[zone_num1].surf_dens*ring_area;
		    ring_area=(dust_budget[zone_num2].rad+0.125)*(dust_budget[zone_num2].rad+0.125)-x*x;
		    ring_area=ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		    budget+=dust_budget[zone_num2].surf_dens*ring_area;
	    }
	    time_tot=time_tot+dt;
            fprintf(fp_vr, "%f\t%f\n",x,vr0);
            //fprintf(fp_drt, "%f\t%f\n",x,-1.0*sum1*AU_km/TUNIT);
            fprintf(fp_size,"%f\t%f\n",x,a_pb1);
            if (a_pb1>9.0/4.0*mean_path(x)) {
              //  l++;
              //  if(l==1) x_stop=x;
              //  coag_eff=coag_eff*exp((x-x_stop)/0.001);
                coag_eff=0.000;
                
            }
            //vol_plus=-1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*step*(k1+2*k2+2*k3+k4)/6.0*AU_km*100000.0;
            vol_plus=1.0*M_PI*a_pb1*a_pb1*sqrt(vr0*vr0+0.25*tau*vr0*tau*vr0)*dt*TUNIT*AU_km;
	   // if (0) vol_plus=0.0;
            a_pb2=pow(((vol_plus*coag_eff*density(x)/rho_peb+4.0/3.0*M_PI*a_pb1*a_pb1*a_pb1)*3.0/4.0/M_PI),0.33333333333333333);
	    massplus=(pow(a_pb2/pp->size[0],3)-pow(a_pb1/pp->size[0],3))*pp->mass[0];
	    
	if(zone_num1==zone_num2){
		ring_area0=(dust_budget[zone_num1].rad+0.125)*(dust_budget[zone_num1].rad+0.125)-(dust_budget[zone_num1].rad-0.125)*(dust_budget[zone_num1].rad-0.125);
		ring_area0=ring_area0*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		ring_mass=dust_budget[zone_num1].surf_dens*ring_area0;
		if(massplus<=budget){
			ring_mass=ring_mass-massplus;
			dust_budget[zone_num1].surf_dens=ring_mass/ring_area0;
		}
		else{
			ring_mass=ring_mass-budget;
			dust_budget[zone_num1].surf_dens=ring_mass/ring_area0;
			a_pb2=pow(budget/massplus,0.333333333333)*a_pb2;
		}
	}
	else{
		ring_area_all=(x+vr0*dt/AU_km*TUNIT)*(x+vr0*dt/AU_km*TUNIT)-x*x;
		ring_area_all=ring_area_all*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		for(i=zone_num1-1;i>zone_num2;i--){
			ring_area=(dust_budget[i].rad+0.125)*(dust_budget[i].rad+0.125)-(dust_budget[i].rad-0.125)*(dust_budget[i].rad-0.125);
	        	ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
			ring_mass=ring_area*dust_budget[i].surf_dens;
			if(massplus<=budget) {
				ring_mass=ring_mass-massplus*ring_area/ring_area_all;
				if (ring_mass<0.0) ring_mass=0.0;
			}
			else ring_mass=0.0;
			dust_budget[i].surf_dens=ring_mass/ring_area;
		}
		ring_area=(x+vr0*dt/AU_km*TUNIT)*(x+vr0*dt/AU_km*TUNIT)-(dust_budget[zone_num1].rad-0.125)*(dust_budget[zone_num1].rad-0.125);
		ring_area=ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		ring_area0=(dust_budget[zone_num1].rad+0.125)*(dust_budget[zone_num1].rad+0.125)-(dust_budget[zone_num1].rad-0.125)*(dust_budget[zone_num1].rad-0.125);
		ring_area0=ring_area0*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
                ring_mass=dust_budget[zone_num1].surf_dens*ring_area0;
		if(massplus<=budget){
			ring_mass=ring_mass-massplus*ring_area/ring_area_all;
			if (ring_mass<0.0) ring_mass=0.0;
	    	        }
		else {
			ring_mass=ring_mass-budget*ring_area/ring_area_all;                                             if (ring_mass<0.0) ring_mass=0.0;
		}
                dust_budget[zone_num1].surf_dens=ring_mass/ring_area;

		ring_area=(dust_budget[zone_num2].rad+0.125)*(dust_budget[zone_num2].rad+0.125)-x*x;
		ring_area=ring_area*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
		ring_area0=(dust_budget[zone_num2].rad+0.125)*(dust_budget[zone_num2].rad+0.125)-(dust_budget[zone_num2].rad-0.125)*(dust_budget[zone_num2].rad-0.125);
		ring_area0=ring_area0*(AU_km*100000.0)*(AU_km*100000.0)*M_PI;
                ring_mass=dust_budget[zone_num2].surf_dens*ring_area0;
		if(massplus<=budget){
			ring_mass=ring_mass-massplus*ring_area/ring_area_all;
			if (ring_mass<0.0) ring_mass=0.0;
	    	        }
		else {
			ring_mass=ring_mass-budget*ring_area/ring_area_all;                                             if (ring_mass<0.0) ring_mass=0.0;
		}
		dust_budget[zone_num2].surf_dens=ring_mass/ring_area;
		if(massplus>budget) 	a_pb2=pow(budget/massplus,0.333333333333)*a_pb2;
	}
		
	    
            pp->vr[ll]=vr0;
	    pp->size[ll+1]=a_pb2;
	    pp->rad[ll+1]=x;
	    pp->time[ll+1]=time_tot;
	    pp->mass[ll+1]=pow(pp->size[ll+1]/pp->size[0],3)*pp->mass[0];
//            fprintf(fp_peb,"%2.12g\t%2.12g\t%2.12g\t%2.12g\n",pp->rad[ll],pp->time[ll],pp->size[ll],pp->mass[ll]);
	    ll++;
    }
    
//    printf("%0.10f\n", sum1);
//    printf("%0.10f\n", tau_fric0(1.0,a_pb1));
//    printf("%0.10f\n", v_K(0.01516));
//    printf("%0.10f\n", yeta(0.01516));
    
    
    
    
}
       pp->max_step=ll-1;
       //printf("MAXSTEP %f\t%f\n",pp->mass[ll],pp->mass[ll-1]);

    fclose(fp_vr);
    fclose(fp_peb);
    fclose(fp_size);
    
    
    sum1=sum1*AU_km/TUNIT;

return 0;
}*/
