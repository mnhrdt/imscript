/*  optmiz.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
#include <stdlib.h>
#define Abs(x) ( ((x)<0.)?-(x):(x) )
static double fev(double *x,double *py,double *ps,
	double c,double (*func)());
int optmiz(double *x,int n,double (*func)(),double de,double test,int max)
{ double fs,fp,fa,fb,fc,s,sa,sb,sc; int k,m;
  double *pd,*ps,*py,*pg,*ph,*p,*q,*r;
  pd=(double *)calloc(n*(n+4),sizeof(double));
  ps=pd+n; py=ps+n; pg=py+n; ph=pg+n;
  for(p=ph,q=ph+n*n; p<q ;p+=n+1) *p=1.;
  for(p=x,q=pg,fb=(*func)(x); q<ph ;){
    *p+=de; *q++ =((*func)(x)-fb)/de; *p++ -=de;}
  for(m=0; m<max ;++m){
    for(p=ps,r=ph; p<py ;++p){ *p=0.;
      for(q=pg; q<ph ;) *p-= *r++ * *q++; }
    fp=fa=fb; sa=sb=0.; sc=1.;
    for(;;){ if((fc=fev(x,py,ps,sc,func))>fb) break;
      fa=fb; sa=sb; fb=fc; sb=sc; sc*=2.; }
    if(sc==1.){ sb=.5;
      for(;;){ if((fb=fev(x,py,ps,sb,func))<fa||sb<1.e-3) break;
         fc=fb; sc=sb; sb/=2.;} }
    for(k=0; k<3 ;++k){ s=(fc-fa)/(sc-sa);
      if((fs=(s-(fb-fa)/(sb-sa))/(sc-sb))<0.) break;
      if((s=(sa+sc-s/fs)/2.)==sb) s-=(sb-sa)/5.;
      fs=fev(x,py,ps,s,func);
      if(fs<fb){ if(s<sb){ sc=sb; fc=fb;} else{ sa=sb; fa=fb;}
         sb=s; fb=fs; }
      else{ if(s<sb){ sa=s; fa=fs;} else{ sc=s; fc=fs;} }
     }
    for(p=x,r=ps; r<py ;){ *r*=sb; *p++ += *r++;}
    if(Abs(fp-fb)<test){ free(pd); return (m+1);}
    for(p=x,q=pg,r=pd; q<ph ;){
      *p+=de; fa=((*func)(x)-fb)/de; *p++ -=de;
      *r++ =fa- *q; *q++ =fa; }
    for(p=py,r=ph; p<pg ;++p){ *p=0.;
      for(q=pd; q<ps ;) *p+= *r++ * *q++; }
    for(p=py,q=ps,r=pd,sa=sb=0.; p<pg ;){
      sa+= *r* *p++; sb+= *q++ * *r++; }
    sa=1.+sa/sb;
    for(p=ps,q=py,r=ph; p<py ;++p,++q){
      for(k=0; k<n ;++k)
        *r++ +=(*(ps+k)* *p*sa- *(py+k)* *p- *(ps+k)* *q)/sb;
       }
   }
  free(pd); return 0;
}
static double fev(double *x,double *py,double *ps,double c,double (*func)())
{ double *p,*q,*r;
  for(p=x,q=py,r=ps; r<py ;) *q++ = *p++ +c* *r++;
  return (*func)(py);
}
