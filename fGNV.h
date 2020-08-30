//Не бейте за этот код. Было поставлено ТЗ использовать данный файл в текущем программном проекте.
#pragma once
#include <cmath>       // Математические функции
double fGNV(short int in,short int iv,double zi,short int& error)
{     // in - первый индекс, iv - второй индекс, zi - действительный аргумент,
    // error - код возврата: 0 - правильный счёт; 1 - недопустимые сочетания аргумента и индексов,
    //                                                не считает, возврат 0.
    short int v,n,nv,k,k1,k2,k3,m;
    double yk,y1,y2,ya,yb,z,z1,z1v,z2,z3,e,ek,en,env,fn,fne,fv,fnv,fve,snv,sos,s1,s2, \
       a,alc,s,d,cc,cn,ck,ck1,ck2;
    yk=0.;         error=0;
    n=in;          v=iv;            z=zi;
    e=1.0e+00;     ek=-1.0e+00;     sos=0.1e-06;
    s1=0.;         s2=0.;           nv=n+v;
    if (nv<0) {n=0; v=1;}
    if (z<0)  {error=1; goto m52;}
    if (z==0.) {if (v<=0) {error=1; goto m52;} else goto m18;}
    if (z<5.5) goto m6;
    if (z-nv>2) goto m35;
    m6:  z1=z/2.;
    z2=z1*z1;
    alc=-(log(z1)+0.577215664901532860606512);
    if (nv>0) goto m9;
    env=e;   fnv=e;   snv=0.;   goto m15;
    m9:  k=1;
    env=e;
    m10:   env=env*ek;
    k=k+1;
    if (k-nv<=0) goto m10;
    k=1;
    fnv=e;
    m12: fnv=fnv*k;
    k=k+1;
    if (k-nv<=0) goto m12;
    k=1;
    snv=0.;
    m14: snv=snv+e/k;
    k=k+1;
    if (k-nv<=0) goto m14;
    m15: ;
    z1v=pow(z1,(2*v));
    a=alc+0.5*snv;
    y1=z1v/fnv;
    k=1;     s=0.;
    m17: d=s;
    s=s+y1*(a+0.5*(s1+s2));
    s1=s1+e/k;
    s2=s2+e/(k+nv);
    y1=y1*z2/(k*(k+nv));
    k=k+1;
    cc=abs(d);    cn=abs(s);
    if (abs(cc-cn)-cc*sos>0.) goto m17;
    ya=env*s;
    m18: ;
    if (v>0) goto m20;
    yb=0.;   goto m34;
    m20: if (n>0) goto m22;
    en=e;     fn=e;     goto m26;
    m22: k=1;
    en=e;
    m23:  en=en*ek;
    k=k+1;
    if (k-n<=0) goto m23;
    k=1;
    fn=e;
    m25:  fn=fn*k;
    k=k+1;
    if (k-n<=0) goto m25;
    m26: ;
    if (v-1<=0) {fve=e; goto m30;}
    k=1;   fve=e;
    m29:  fve=fve*k;
    k=k+1;
    if (k-v<0) goto m29;
    m30:  y2=fve/fn;
    k=0;     s=0.;
    m33:  s=s+y2;
    k=k+1;
    if (k-v+1<=0) {y2=y2*ek*z2/((v-k)*(n+k));   goto m33;}
    yb=en*s;
    m34:  yk=ya+0.5*yb;   goto m52;
    // Расчёт по асимптотической формуле.
    m35:  z1=2./z;   z2=z1*z1;   z3=0.125/z;   k1=nv+5;   m=nv*nv*4;
    y1=1.;   s=1.;   k=1;
    m36:  k2=2*k-1;
    k3=k2*k2;
    y1=y1*(m-k3)*z3/k;   s=s+y1;   k=k+1;
    if (k-k1<=0) goto m36;
    ck1=pow(z/2.,v-n);
    ck2=1.2533141/sqrt(z);
    ck=ck1*ck2/exp(z);
    ya=ck*s;
    if (n<=0) {yb=0.;  goto m51;}
    k=1;   en=e;
    m40:  en=en*ek;
    k=k+1;
    if (k-n<=0) goto m40;
    k=1;   fv=e;
    m42:  fv=fv*k;
    k=k+1;
    if (k-n<=0) goto m42;
    if (n-1<=0) {fne=e;  goto m47;}
    k=1;   fne=e;
    m46:  fne=fne*k;
    k=k+1;
    if (k-n<=0) goto m46;
    m47:  y2=ek*fv*z2/fne;
    k=1;   s=0.;
    m50:  s=s+y2;
    if (k-n<0) {y2=y2*ek*(v+k)*(n-k)*z2;  k=k+1;  goto m50;}
    yb=en*s;
    m51:  yk=ya-0.5*yb;
    m52:  ;
    return yk;
}