#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
#include<math.h>
#include<stdio.h>
void frag(){
int i,j,jl,js;
double size_lim,frac,AREA,dr,delta_mass;
dr=size_ring;
for(i=ring_num-1;i>=0;i--){
	dr=peb_map[i].dr;
	for(j=peb_size_num-1;j>=0;j--){
		size_lim=2.25*mean_path(peb_map[i].rad);		
		if(peb_map[i].size[j] > 0.8*size_lim){
			jl=j;
			js=floor(log10(peb_map[i].size[jl]/0.1/2.0)/size_step);
//			js=jl;
			if(js <0) js=0;
			frac=peb_map[i].size[jl]/size_lim;
			frac=frac*0.9;
			if (frac >1.0) frac=1.0;
			//frac=0.0;
			delta_mass=peb_map[i].mass_out[j]*frac;
			peb_map[i].mass_out[jl]-=delta_mass;
			peb_map[i].mass_out[js]+=delta_mass;
	}
}
}

for(i=0;i<ring_num;i++){
	AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
	for(j=0;j<peb_size_num;j++){
		peb_map[i].surf_dens[j]=(peb_map[i].mass_out[j]+peb_map[i].mass_in[j])/AREA;
	}
}
}

double coag_kernel(double a_pb1,double a_pb2,double delta_v,double n_density1, double n_density2,double hei1, double hei2, int i, int j, int jj){
	double cross_section, path,n_density;
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*n_density1*n_density2/sqrt(2*M_PI*(hei1*hei1+hei2*hei2));
//	if(delta_v>v_coag_max+v_tran_width) path=0.0;
//	else if(delta_v>v_coag_max && delta_v<=v_coag_max+v_tran_width) path\
		=(v_coag_max+v_tran_width-delta_v)/v_tran_width*path;
	//if(0 || path/n_density2>1.0)	printf("dv=%e path=%e N1=%e  N2=%e  sig1=%e hei=%e %d %d %d %e %e\n",\
			delta_v/100,path/n_density2,n_density1,n_density2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj,peb_map[i].size[j],peb_map[i].size[jj]);
      //printf("PATH=%e\t",path);
	return 1.0*path;
}


double coag_kernel_sweep(double a_pb1,double a_pb2,double delta_v,double n_density1, double n_density2,double hei1, double hei2, int i, int j, int jj){
	double cross_section, path,n_density;
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*n_density1*n_density2/sqrt(2*M_PI*(hei1*hei1+hei2*hei2));
	if(path>1.0)	printf("path=%e N1=%e  N2=%e  sig1=%e hei=%e %d %d %d\n",\
			path,n_density1,n_density2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj);
	return 1.0*path;
}


double coag_kernel2(double a_pb1,double a_pb2,double delta_v,double rho1, double rho2, int i, int j, int jj){
	double cross_section, path,n_density,n_density1,n_density2;
	n_density2=rho2/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0)*0.85;
	n_density1=rho1/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0)*0.85;//disk within -h to h
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	//cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*rho2;
        if(delta_v>v_coag_max+v_tran_width) path=0.0;
	else if(delta_v>v_coag_max && delta_v<=v_tran_width) path\
	       	=(v_coag_max+v_tran_width-delta_v)/v_tran_width*path;
	if(path>1.0)	printf("path=%e N1=%e rho1=%e N2=%e rho2=%e sig1=%e hei=%e %d %d %d\n",\
			path,n_density1,rho1,n_density2,rho2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj);
	return 1.0*path;
}

double frag_size_func(double size){// size in cm
  //  return exp(-1.0*size_slope*peb_map[i].size_med[j]);
    //return pow(size,-3.0*11.0/6.0/size_slope_pow+5.0);
    return pow(size,frag_slope);
}

void coagulation(double dt0, double tot_time){
        int i,j,jj,jjj=0,j_new,j_news,j_newl,N_diff=100000,sweep=0,i_inter_min;
	double a_pb1,a_pb2,a_pb3,v1,v2,delta_v,dr,dt1,sub_t;
	double frac_s,frac_s2,frac_m,frac_m1,frac_m2,mass0,n1,n2,n0,n_delta,h1,h2,m_pb1,m_pb2,m_pb3,AREA;
	double frac_m1t,frac_m2t;//temporary mass fraction
	
	if(1 && ((int)tot_time)%100==0){
	for(i=0;i<ring_num;i++){
	      dtco_ring[i]=dt0;
		}
	}
  for(i=0;i<ring_num;i++){
    if (peb_map[i].rad>=interaction_min_radius)
    {
      i_inter_min=i;
      break;
    }
  }
  //i_inter_min=1;
	for (i=ring_num-1;i>=i_inter_min;i--){
		dr=peb_map[i].dr;
		AREA=peb_map[i].AREA;
		frac_s=0.0;
		frac_m=0.0;
		frac_m1=0.0;
		frac_m2=0.0;
		frac_s2=0.0;
		

	if(1 && ((int)tot_time)%100!=0){
			 dt1=2.0*dtco_ring[i];
	}
	else {
		dt1=2.0*dt0;
	}
		
		do{
		dt1=dt1/2.0;
		frac_s=0.0;
		frac_m=0.0;
		frac_m1=0.0;
		frac_m2=0.0;
		for(j=peb_size_num-1;j>=1;j--){
			if(peb_map[i].surf_dens[j] < 1e-5) continue;
			for(jj=j-1;jj>=0;jj--){
			if(peb_map[i].surf_dens[jj] < 1e-5) continue;
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];
			n1=peb_map[i].surf_dens[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			n2=peb_map[i].surf_dens[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
      v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			if(0||jj>=(((j-N_diff)>0)?(j-N_diff):0)){
			n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
      if(1 && FRAG_SW){
        if(delta_v>v_frag_max0) n_delta=0.0;
        else if(delta_v>=v_frag_max1 && delta_v<=v_frag_max0) {
        n_delta=(v_frag_max0-delta_v)/(v_frag_max0-v_frag_max1)*n_delta;
        }
        else if(delta_v<v_frag_min1 && delta_v>=v_frag_min0) {
        n_delta=(delta_v-v_frag_min0)/(v_frag_min1-v_frag_min0)*n_delta;
        }
        else if( delta_v<v_frag_min0 && delta_v>=v_coag_max+v_tran_width) n_delta=0.0;
      }
			if(!FRAG_SW && delta_v>v_coag_max+v_tran_width) n_delta=0.0;
      else if(delta_v>v_coag_max && delta_v<=v_coag_max+v_tran_width) n_delta\
                =(v_coag_max+v_tran_width-delta_v)/v_tran_width*n_delta;




//			if(n1<n2) frac_m2=n_delta/n1;//(2*M_PI*h1*h2*n2/sqrt(2*M_PI*(h1*h1+h2*h2)));
//			else frac_m2=n_delta/n2;//(2*M_PI*h1*h2*n1/sqrt(2*M_PI*(h1*h1+h2*h2)));

			frac_m1t=n_delta*m_pb2/n1/m_pb1;
			frac_m2t=n_delta/n2;
			if(frac_m1t>frac_m1) frac_m1=frac_m1t;
			if(frac_m2t>frac_m2) frac_m2=frac_m2t;

			}
			else if(sweep){
				//frac_m1=1.0/m_pb1*dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],i,j,jj);
				frac_m1=1.0/n1/m_pb1*m_pb2*dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);	
				//frac_m1=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
				frac_m2=frac_m1*n1*m_pb1/peb_map[i].mass_out[jj];
				frac_s2=pow(frac_m1,1.0/3.0);
				if(frac_m2<frac_m1) frac_m2=frac_m1;
					
				
	//			if(frac_s2>0.1) printf("AHAGAGA Frac_S=%e i=%d j=%d jj=%d\n",frac_s2,i,j,jj);
			}



			//frac_s2=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/n2*h2;
			//else frac_s2=n_delta/n1*h1;
			
			if(frac_m2>frac_m) frac_m=frac_m2;

			if(frac_s2>frac_s){
			frac_s=frac_s2;
			jjj=j;
			}
			
			}
			}
			//printf("frac_m=%e\t",frac_m);	
			}while(frac_m1>0.02 ||  frac_m2>0.02);
		dtco_ring[i]=dt1;
	if( 0 && dt1 < dt0) printf("ring=%d\tdt=%g %f max_frac_s=%f j=%d\n",i,dt1,dt0,frac_s,jjj);
	sub_t=0.0;
	while(sub_t<dt0){
		for(j=peb_size_num-1;j>=1;j--){// original sequence
		//for(j=0;j<peb_size_num;j++){// test sequence
      if(peb_map[i].surf_dens[j] < 1e-5) continue;
      for(jj=j-1;jj>=0;jj--){// original sequence
      //for(jj=j+1;jj<peb_size_num;jj++){// test sequence

			
			
				
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];

			
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			n1=peb_map[i].surf_dens[j]/m_pb1;
			n2=peb_map[i].surf_dens[jj]/m_pb2;

			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
                        v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			
			if(0||jj>=(((j-N_diff)>0)?(j-N_diff):0)){
			double frac_ss,frac_sm,frac_sl;
			double a_pb1s,a_pb1l,a_pb2s,a_pb2l,a_pb3s,a_pb3l;	
			int j_bet;	
				a_pb1s=peb_map[i].size[j];
				a_pb1l=peb_map[i].size[j+1];
				a_pb2s=peb_map[i].size[jj];
				a_pb2l=peb_map[i].size[jj+1];
				a_pb3s=pow(a_pb1s*a_pb1s*a_pb1s+a_pb2s*a_pb2s*a_pb2s,1.0/3.0);
				a_pb3l=pow(a_pb1l*a_pb1l*a_pb1l+a_pb2l*a_pb2l*a_pb2l,1.0/3.0);
				a_pb3=pow(a_pb1*a_pb1*a_pb1+a_pb2*a_pb2*a_pb2,1.0/3.0);
				
				j_new=floor(log10(a_pb3/size_min)/size_step);	
				j_newl=floor(log10(a_pb3l/size_min)/size_step);
				j_news=floor(log10(a_pb3s/size_min)/size_step);
					
				if(j_new>=peb_size_num) j_new=peb_size_num-1;
				if(j_newl>=peb_size_num) j_newl=peb_size_num-1;
				if(j_news>=peb_size_num) j_news=peb_size_num-1;
				
				frac_ss=(peb_map[i].size[j_news+1]-a_pb3s)/(a_pb3l-a_pb3s);
				//frac_sm=(peb_map[i].size[j_new]-peb_map[i].size[j_news+1])\
					/(a_pb3l-a_pb3s);  
	
				frac_sl=(a_pb3l-peb_map[i].size[j_newl])\
					/(a_pb3l-a_pb3s);   
				frac_sm=1-frac_ss-frac_sl;

				//m_pb3=rho_peb*4.0*M_PI*a_pb3*a_pb3*a_pb3/3.0;
				m_pb3=m_pb1+m_pb2;
				n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
        if(1 && FRAG_SW){// && delta_v>v_coag_max+v_tran_width)
          int j_frag=0;
          double frag_mass_norm=0.0,mass_collide=0.0,mass_frag=0.0;
          if(delta_v>=v_frag_max0) n_delta=0.0;
          else if(delta_v>=v_frag_max1 && delta_v<v_frag_max0) {
           n_delta=(v_frag_max0-delta_v)/(v_frag_max0-v_frag_max1)*n_delta;
          }
          else if(delta_v<v_frag_min1 && delta_v>v_frag_min0) {
           n_delta=(delta_v-v_frag_min0)/(v_frag_min1-v_frag_min0)*n_delta;
           }
          else if( delta_v<=v_frag_min0) n_delta=0.0;
//n_delta=0.0;
          /*
          for( j_frag=0;j_frag<jj;j_frag++){        
          //  frag_mass_norm+=p_size_func(peb_map[i].size_med[j_frag]);
            frag_mass_norm+=p_size_func(size_min*pow(10,j_frag*size_step));
          }
          for( j_frag=0;j_frag<1;j_frag++){
            if(0 && j_frag>=0){
            //peb_map[i].mass_in[j_frag]+=m_pb3*n_delta*AREA*frag_size_func(size_min*pow(10,j_frag*size_step))/frag_mass_norm;
            //peb_map[i].mass_in[j_frag]+=m_pb3*n_delta*AREA*frag_size_func(size_min*pow(10,j_frag*size_step))/frag_mass_norm;
            peb_map[i].mass_in[j_frag]+=0.0*m_pb3*n_delta*AREA/(jj);
            mass_collide+=m_pb3*n_delta*AREA/(jj);
            }
            else{
              //dust_budget[i].mass_out+=m_pb3*n_delta*AREA*frag_size_func(peb_map[i].size_med[j_frag])/frag_mass_norm;
              dust_budget[i].mass_out+=m_pb3*n_delta*AREA;
              dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
              dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
            }
          }
         */

              dust_budget[i].mass_out+=(m_pb1+m_pb2)*n_delta*AREA;
              dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
              dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
            
            peb_map[i].mass_in[jj]-=m_pb2*n_delta*AREA;
            peb_map[i].mass_in[j]-=m_pb1*n_delta*AREA;
            mass_frag=m_pb2*n_delta*AREA+m_pb1*n_delta*AREA;
            //if(n_delta>0.0) printf("MASS_COL=%e\tMASS_FRAG=%e\t",mass_collide,mass_frag);

        }

        n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
      //  n_delta=0.0;
        if(delta_v>v_coag_max+v_tran_width) n_delta=0.0;
        else if(delta_v>v_coag_max && delta_v<=v_coag_max+v_tran_width) n_delta\
          =(v_coag_max+v_tran_width-delta_v)/v_tran_width*n_delta;


				if(j_newl-j_news-1>0){	
				for(j_bet=j_news+1;j_bet<j_newl;j_bet++){
					peb_map[i].mass_in[j_bet]+=frac_sm*m_pb3*n_delta*AREA/(j_newl-j_news-1);
				}
				}
				else	peb_map[i].mass_in[j_new]+=frac_sm*m_pb3*n_delta*AREA;
				
				peb_map[i].mass_in[j_news]+=frac_ss*m_pb3*n_delta*AREA;
				peb_map[i].mass_in[j_newl]+=frac_sl*m_pb3*n_delta*AREA;

  //test for frag in coag mode
 //       dust_budget[i].mass_out+=(m_pb1+m_pb2)*n_delta*AREA;
  //      dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
   //     dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
                        	peb_map[i].mass_in[jj]-=m_pb2*n_delta*AREA;
                        	peb_map[i].mass_in[j]-=m_pb1*n_delta*AREA;			
			}
			else if(sweep){	
				j_new=j+1;
				if(j_new>peb_size_num-1) j_new=j;
				//frac_m1=1.0/m_pb1*dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],i,j,jj);
				frac_m1=1.0/n1/m_pb1*m_pb2*dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
				
                                frac_m2=frac_m1*n1*m_pb1/peb_map[i].mass_out[jj];
                                frac_s=pow(frac_m1,1.0/3.0);
				//frac_s=pow(m_pb2*n_delta/n1/m_pb1+1,1.0/3.0)-1.0;
				peb_map[i].mass_in[jj]-=frac_m1*n1*AREA;
				peb_map[i].mass_in[j]-=frac_s*peb_map[i].mass_out[j];
				peb_map[i].mass_in[j_new]+=frac_m1*n1*AREA+frac_s*peb_map[i].mass_out[j];
				
			}

			//frac_s=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/(2*M_PI*h1*h2*n2/sqrt(2*M_PI(h1*h1+h2*h2)));
			//else frac_s2=n_delta/(2*M_PI*h1*h2*n1/sqrt(2*M_PI(h1*h1+h2*h2)));
			
			
			//		if(frac_s> 0.8) printf("rho=%e\t%e\t%f\t%f\t%f\t%f\t%f\n",peb_map[i].rho[j],peb_map[i].rho[jj],frac_s,peb_map[i].rad_med,peb_map[i].size_med[j],peb_map[i].size_med[jj],delta_v);
//			peb_map[i].mass_in[j_new]+=mass0*frac_s*(pow(a_pb2/a_pb1,3)+1.0);
//			peb_map[i].mass_out[jj]-=mass0*frac_s*pow(a_pb2/a_pb1,3);
//			peb_map[i].mass_out[j]-=mass0*frac_s;
			
			}
			
		}
		for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        	}
		
	sub_t+=dt1;

	}
	}

        for(i=ring_num-1;i>=i_inter_min;i--){
        for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        }
        }


}
void fragmentation(double dt0, double tot_time){
        int i,j,jj,jjj=0,j_new,j_news,j_newl,N_diff=100000,sweep=0,i_inter_min;
	double a_pb1,a_pb2,a_pb3,v1,v2,delta_v,dr,dt1,sub_t;
	double frac_s,frac_s2,frac_m,frac_m1,frac_m2,mass0,n1,n2,n0,n_delta,h1,h2,m_pb1,m_pb2,m_pb3,AREA;
	double frac_m1t,frac_m2t;//temporary mass fraction
	
	if(1 && ((int)tot_time)%100==0){
	for(i=0;i<ring_num;i++){
	      dtco_ring[i]=dt0;
		}
	}
  for(i=0;i<ring_num;i++){
    if (peb_map[i].rad>=interaction_min_radius)
    {
      i_inter_min=i;
      break;
    }
  }
  //i_inter_min=1;
	for (i=ring_num-1;i>=i_inter_min;i--){
		dr=peb_map[i].dr;
		AREA=peb_map[i].AREA;
		frac_s=0.0;
		frac_m=0.0;
		frac_m1=0.0;
		frac_m2=0.0;
		frac_s2=0.0;
		

	if(1 && ((int)tot_time)%100!=0){
			 dt1=2.0*dtco_ring[i];
	}
	else {
		dt1=2.0*dt0;
	}
		
		do{
		dt1=dt1/2.0;
		frac_s=0.0;
		frac_m=0.0;
		frac_m1=0.0;
		frac_m2=0.0;
		for(j=peb_size_num-1;j>=1;j--){
			if(peb_map[i].surf_dens[j] < 1e-5) continue;
			for(jj=j-1;jj>=0;jj--){
			if(peb_map[i].surf_dens[jj] < 1e-5) continue;
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];
			n1=peb_map[i].surf_dens[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			n2=peb_map[i].surf_dens[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
      v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			if(0||jj>=(((j-N_diff)>0)?(j-N_diff):0)){
			n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
      if(1 && FRAG_SW){
        if(delta_v>v_frag_max0) n_delta=0.0;
        else if(delta_v>=v_frag_max1 && delta_v<=v_frag_max0) {
        n_delta=(v_frag_max0-delta_v)/(v_frag_max0-v_frag_max1)*n_delta;
        }
        else if(delta_v<v_frag_min1 && delta_v>=v_frag_min0) {
        n_delta=(delta_v-v_frag_min0)/(v_frag_min1-v_frag_min0)*n_delta;
        }
        else if( delta_v<v_frag_min0 && delta_v>=v_coag_max+v_tran_width) n_delta=0.0;
      }
			if(!FRAG_SW && delta_v>v_coag_max+v_tran_width) n_delta=0.0;
      else if(delta_v>v_coag_max && delta_v<=v_coag_max+v_tran_width) n_delta\
                =(v_coag_max+v_tran_width-delta_v)/v_tran_width*n_delta;




//			if(n1<n2) frac_m2=n_delta/n1;//(2*M_PI*h1*h2*n2/sqrt(2*M_PI*(h1*h1+h2*h2)));
//			else frac_m2=n_delta/n2;//(2*M_PI*h1*h2*n1/sqrt(2*M_PI*(h1*h1+h2*h2)));

			frac_m1t=n_delta*m_pb2/n1/m_pb1;
			frac_m2t=n_delta/n2;
			if(frac_m1t>frac_m1) frac_m1=frac_m1t;
			if(frac_m2t>frac_m2) frac_m2=frac_m2t;

			}
			else if(sweep){
				//frac_m1=1.0/m_pb1*dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],i,j,jj);
				frac_m1=1.0/n1/m_pb1*m_pb2*dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);	
				//frac_m1=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
				frac_m2=frac_m1*n1*m_pb1/peb_map[i].mass_out[jj];
				frac_s2=pow(frac_m1,1.0/3.0);
				if(frac_m2<frac_m1) frac_m2=frac_m1;
					
				
	//			if(frac_s2>0.1) printf("AHAGAGA Frac_S=%e i=%d j=%d jj=%d\n",frac_s2,i,j,jj);
			}



			//frac_s2=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/n2*h2;
			//else frac_s2=n_delta/n1*h1;
			
			if(frac_m2>frac_m) frac_m=frac_m2;

			if(frac_s2>frac_s){
			frac_s=frac_s2;
			jjj=j;
			}
			
			}
			}
			//printf("frac_m=%e\t",frac_m);	
			}while(frac_m1>0.02 ||  frac_m2>0.02);
		dtco_ring[i]=dt1;
	if( 0 && dt1 < dt0) printf("ring=%d\tdt=%g %f max_frac_s=%f j=%d\n",i,dt1,dt0,frac_s,jjj);
	sub_t=0.0;
	while(sub_t<dt0){
		for(j=peb_size_num-1;j>=1;j--){// original sequence
		//for(j=0;j<peb_size_num;j++){// test sequence
      if(peb_map[i].surf_dens[j] < 1e-5) continue;
      for(jj=j-1;jj>=0;jj--){// original sequence
      //for(jj=j+1;jj<peb_size_num;jj++){// test sequence

			
			
				
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];

			
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			n1=peb_map[i].surf_dens[j]/m_pb1;
			n2=peb_map[i].surf_dens[jj]/m_pb2;

			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
                        v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			
			if(0||jj>=(((j-N_diff)>0)?(j-N_diff):0)){
			double frac_ss,frac_sm,frac_sl;
			double a_pb1s,a_pb1l,a_pb2s,a_pb2l,a_pb3s,a_pb3l;	
			int j_bet;	
				a_pb1s=peb_map[i].size[j];
				a_pb1l=peb_map[i].size[j+1];
				a_pb2s=peb_map[i].size[jj];
				a_pb2l=peb_map[i].size[jj+1];
				a_pb3s=pow(a_pb1s*a_pb1s*a_pb1s+a_pb2s*a_pb2s*a_pb2s,1.0/3.0);
				a_pb3l=pow(a_pb1l*a_pb1l*a_pb1l+a_pb2l*a_pb2l*a_pb2l,1.0/3.0);
				a_pb3=pow(a_pb1*a_pb1*a_pb1+a_pb2*a_pb2*a_pb2,1.0/3.0);
				
				j_new=floor(log10(a_pb3/size_min)/size_step);	
				j_newl=floor(log10(a_pb3l/size_min)/size_step);
				j_news=floor(log10(a_pb3s/size_min)/size_step);
					
				if(j_new>=peb_size_num) j_new=peb_size_num-1;
				if(j_newl>=peb_size_num) j_newl=peb_size_num-1;
				if(j_news>=peb_size_num) j_news=peb_size_num-1;
				
				frac_ss=(peb_map[i].size[j_news+1]-a_pb3s)/(a_pb3l-a_pb3s);
				//frac_sm=(peb_map[i].size[j_new]-peb_map[i].size[j_news+1])\
					/(a_pb3l-a_pb3s);  
	
				frac_sl=(a_pb3l-peb_map[i].size[j_newl])\
					/(a_pb3l-a_pb3s);   
				frac_sm=1-frac_ss-frac_sl;

				//m_pb3=rho_peb*4.0*M_PI*a_pb3*a_pb3*a_pb3/3.0;
				m_pb3=m_pb1+m_pb2;
				n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
        if(1 && FRAG_SW){// && delta_v>v_coag_max+v_tran_width)
          int j_frag=0;
          double frag_mass_norm=0.0,mass_collide=0.0,mass_frag=0.0;
          if(delta_v>=v_frag_max0) n_delta=0.0;
          else if(delta_v>=v_frag_max1 && delta_v<v_frag_max0) {
           n_delta=(v_frag_max0-delta_v)/(v_frag_max0-v_frag_max1)*n_delta;
          }
          else if(delta_v<v_frag_min1 && delta_v>v_frag_min0) {
           n_delta=(delta_v-v_frag_min0)/(v_frag_min1-v_frag_min0)*n_delta;
           }
          else if( delta_v<=v_frag_min0) n_delta=0.0;
//n_delta=0.0;
          /*
          for( j_frag=0;j_frag<jj;j_frag++){        
          //  frag_mass_norm+=p_size_func(peb_map[i].size_med[j_frag]);
            frag_mass_norm+=p_size_func(size_min*pow(10,j_frag*size_step));
          }
          for( j_frag=0;j_frag<1;j_frag++){
            if(0 && j_frag>=0){
            //peb_map[i].mass_in[j_frag]+=m_pb3*n_delta*AREA*frag_size_func(size_min*pow(10,j_frag*size_step))/frag_mass_norm;
            //peb_map[i].mass_in[j_frag]+=m_pb3*n_delta*AREA*frag_size_func(size_min*pow(10,j_frag*size_step))/frag_mass_norm;
            peb_map[i].mass_in[j_frag]+=0.0*m_pb3*n_delta*AREA/(jj);
            mass_collide+=m_pb3*n_delta*AREA/(jj);
            }
            else{
              //dust_budget[i].mass_out+=m_pb3*n_delta*AREA*frag_size_func(peb_map[i].size_med[j_frag])/frag_mass_norm;
              dust_budget[i].mass_out+=m_pb3*n_delta*AREA;
              dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
              dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
            }
          }
         */

              dust_budget[i].mass_out+=(m_pb1+m_pb2)*n_delta*AREA;
              dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
              dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
            
            peb_map[i].mass_in[jj]-=m_pb2*n_delta*AREA;
            peb_map[i].mass_in[j]-=m_pb1*n_delta*AREA;
            mass_frag=m_pb2*n_delta*AREA+m_pb1*n_delta*AREA;
            //if(n_delta>0.0) printf("MASS_COL=%e\tMASS_FRAG=%e\t",mass_collide,mass_frag);

        }

          if(delta_v>=v_frag_max0) n_delta=0.0;
          else if(delta_v>=v_frag_max1 && delta_v<v_frag_max0) {
            n_delta=(v_frag_max0-delta_v)/(v_frag_max0-v_frag_max1)*n_delta;                                           }
          else if(delta_v<v_frag_min1 && delta_v>v_frag_min0) {
            n_delta=(delta_v-v_frag_min0)/(v_frag_min1-v_frag_min0)*n_delta;                                                      }
          else if( delta_v<=v_frag_min0) n_delta=0.0;

/*
				if(j_newl-j_news-1>0){	
				for(j_bet=j_news+1;j_bet<j_newl;j_bet++){
					peb_map[i].mass_in[j_bet]+=frac_sm*m_pb3*n_delta*AREA/(j_newl-j_news-1);
				}
				}
				else	peb_map[i].mass_in[j_new]+=frac_sm*m_pb3*n_delta*AREA;
				
				peb_map[i].mass_in[j_news]+=frac_ss*m_pb3*n_delta*AREA;
				peb_map[i].mass_in[j_newl]+=frac_sl*m_pb3*n_delta*AREA;
*/
  //test for frag in coag mode
        dust_budget[i].mass_out+=(m_pb1+m_pb2)*n_delta*AREA;
        dust_budget[i].surf_dens=dust_budget[i].mass_out/AREA;
        dust_budget[i].rho=dust_budget[i].surf_dens/dust_budget[i].hei;
                        	peb_map[i].mass_in[jj]-=m_pb2*n_delta*AREA;
                        	peb_map[i].mass_in[j]-=m_pb1*n_delta*AREA;			
			}
			else if(sweep){	
				j_new=j+1;
				if(j_new>peb_size_num-1) j_new=j;
				//frac_m1=1.0/m_pb1*dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],i,j,jj);
				frac_m1=1.0/n1/m_pb1*m_pb2*dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
				
                                frac_m2=frac_m1*n1*m_pb1/peb_map[i].mass_out[jj];
                                frac_s=pow(frac_m1,1.0/3.0);
				//frac_s=pow(m_pb2*n_delta/n1/m_pb1+1,1.0/3.0)-1.0;
				peb_map[i].mass_in[jj]-=frac_m1*n1*AREA;
				peb_map[i].mass_in[j]-=frac_s*peb_map[i].mass_out[j];
				peb_map[i].mass_in[j_new]+=frac_m1*n1*AREA+frac_s*peb_map[i].mass_out[j];
				
			}

			//frac_s=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/(2*M_PI*h1*h2*n2/sqrt(2*M_PI(h1*h1+h2*h2)));
			//else frac_s2=n_delta/(2*M_PI*h1*h2*n1/sqrt(2*M_PI(h1*h1+h2*h2)));
			
			
			//		if(frac_s> 0.8) printf("rho=%e\t%e\t%f\t%f\t%f\t%f\t%f\n",peb_map[i].rho[j],peb_map[i].rho[jj],frac_s,peb_map[i].rad_med,peb_map[i].size_med[j],peb_map[i].size_med[jj],delta_v);
//			peb_map[i].mass_in[j_new]+=mass0*frac_s*(pow(a_pb2/a_pb1,3)+1.0);
//			peb_map[i].mass_out[jj]-=mass0*frac_s*pow(a_pb2/a_pb1,3);
//			peb_map[i].mass_out[j]-=mass0*frac_s;
			
			}
			
		}
		for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        	}
		
	sub_t+=dt1;

	}
	}

        for(i=ring_num-1;i>=i_inter_min;i--){
        for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        }
        }


}
