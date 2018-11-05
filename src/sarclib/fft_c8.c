/***************************************************************************
 * 
 *           Module :  fft_c8.c
 *          Program :  sarclib
 *        Created by: Andrew Horne October 1980
 *                    Alan Blake      26/1/1993
 *                    Darren Muff     22/2/1995
 *                    Matt Nottingham 12/09/2005.
 *
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *       C version of Singleton FFT algorithm converted from ADLIB and
 *       FORTRAN code used at RSRE in the 1980's and 1990's. Historical
 *       information attached to explain old-school fft algorithm
 *
 ***************************************************************************
 *
 * REV 0.0                            --- ADLIB --- OCTOBER 1980 ---
 *
 *	Frigged by A.M. Horne to work with a complex array
 *
 *	Call using FFT(SIG(0),SIG(0),...) where the first element is in
 *	SIG(0)
 *
 *        PURPOSE
 *           TO COMPUTE A COMPLEX(MULTIVARIATE)FOURIER TRANSFORM USING A
 *           GENERALIZATION OF THE COOLEY AND TUKEY ALGORITHM,THE INVERSE
 *           TRANSFORMATION CAN ALSO BE DONE
 *
 *        USAGE
 *           CALL FFT(A,B,NTOT,N,NSPAN,ISN)
 *
 *        DESCRIPTION OF PARAMETERS
 *           A      -REAL ARRAY CONTAINING THE REAL COMPONENTS OF THE
 *                   DATA OR THE REAL COMPONENTS OF THE TRANSFORM
 *           B      -REAL ARRAY CONTAINING THE IMAGINARY COMPONENTS OF
 *                   THE DATA OR THE IMAGINARY COMPONENTS OF THE TRANSFORM
 *                   A AND B MAY BE MULTIVARIATE
 *           NTOT   -INTEGER IS THE TOTAL NUMBER OF COMPLEX DATA VALUES
 *           N      -INTEGER IS THE DIMENSION OF THE CURRENT VARIABLE
 *           NSPAN  -INTEGER,NSPAN/N IS THE SPACING OF CONSECTIVE DATA
 *                   VALUES WHILE INDEXING THE CURRENT VARIABLE
 *                   THE MEANING OF NTOT,N AND NSPAN MAY BE UNDERSTOOD
 *                   FROM THE FOLLOWING EXAMPLES WHERE IT IS IMPORTANT TO
 *                   KNOW THAT IN FORTRAN THE ELEMENT A(2,J,...) IS FOUND
 *                   NEXT TO A(1,J,...)
 *                   A  SINGLE-VARIATE TRANSFORM WITH A(N) AND B(N) IS
 *                   COMPUTED BY ONE CALL
 *                             CALL FFT(A,B,N,N,N,ISN)
 *
 *                   A BI-VARIATE TRANSFORM WITH A(N1,N2) AND B(N1,N2) IS
 *                   COMPUTED BY TWO CALLS
 *                             CALL FFT(A,B,N1*N2,N1,N1,ISN)
 *                             CALL FFT(A,B,N1*N2,N2,N1*N2,ISN)
 *
 *                   A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3) AND
 *                   B(N1,N2,N3) IS COMPUTED BY THREE CALLS
 *                             CALL FFT(A,B,N1*N2*N3,N1,N1      ,ISN)
 *                             CALL FFT(A,B,N1*N2*N3,N2,N1*N2   ,ISN)
 *                             CALL FFT(A,B,N1*N2*N3,N3,N1*N2*N3,ISN)
 *
 *           ISN    -INTEGER DEFINES THE TYPE OF TRANSFORM
 *                    ISN=-1 SIGN OF EXPONENT IS NEGATIVE
 *                    ISN= 1 SIGN OF EXPONENT IS POSITIVE
 *
 *        REMARKS
 *           I) THE ROUTINE CAN BE USED WHEN THE DATA IS KEPT IN A SINGLE
 *              COMPLEX ARRAY C,THE PROPER CALL IN THIS CASE IS
 *                             CALL FFT(C,C(2),NTOT,N,NSPAN,2)
 *
 *              THE MAGNITUDE OF ISN IS CHANGED IN TWO
 *          II) IF THE AVAILABLE STORAGE IS INSUFFICIENT ,THE FOLLOWING
 *              MESSAGE IS ISSUED
 *               "ARRAY BOUNDS EXCEEDED WITHIN SUBROUTINE FFT"
 *              IT MEANS THAT
 *                 .THE MAXIMUM PRIMEFACTOR IN N IS GREATER THAN 23
 *                 .THE NUMBER OF PRIMEFACTORS IN N IS GREATER THAN 209
 *         III) THE ROUTINE IS SUBMITTED BY H.OUDSHOORN,UTRECHT(ACCU)
 *
 *        SUBROUTINES AND FUNCTIONS REQUIRED
 *           NONE
 *
 *        METHOD
 *           THE METHOD AND A PROGRAM LISTING ARE DESCRIBED IN THE
 *           FOLLOWING ARTICLE
 *           R.C.SINGLETON,AN ALGORITHM FOR COMPUTING THE MIXED RADIX
 *           FAST FOURIER TRANSFORM,
 *           IEEE TRANS.AUDIO AND ELECTROACOUSTICS VOL AU-17.PP.93-103
 *           JUNE 1969
 *
 *        COMPUTER SYSTEM INFORMATION
 *           TYPE              CYBER 73
 *           SCOPE VERSION     3.4
 *           COMPILER VERSION  FTN 4.0+P357
 *           OPTION USED       OPT=1
 *           DATE              021574
 *
 ***************************************************************************
 *
 *   Copyright (c) 1980 Royal Signals and Radar Establishment. All rights reserved.
 *   Copyright (c) 1993,1996 Defence Research Agency, Malvern. All rights reserved.
 *   Copyright (c) 2013 [dstl]. All rights reserved. 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 ***************************************************************************/


#include "sarclib.h"

#ifndef ROUND
#ifndef TRUNCATE
MUST DEFINE ROUND OR TRUNCATE
#endif
#endif

int fft_c8(SPCmplx *data, int ntot, int n, int nspan, int isn)
{
	float *a, *b;
	int nfac[24],np[210];
    // array storage for maximum prime factor of 61
    //
	float at[62],ck[62],bt[62],sk[62];
	int maxf,maxp,inc,nt,ks,kspan,nn,jc,i,jf;
	int m,k,j,jj,kt,kk,k1,k2,k3,k4,kspnn;
    
	float rad,s72,c72,s120,radf,sd,cd,ak,bk,c1,s1,aj,bj;
	float akp,akm,ajp,ajm,bkp,bkm,bjp,bjm,aa,bb;
	float c2,s2,c3,s3;
    
    // The following two lines stop compiler warnings
    // about these variables being used uninitialised
    // in the code. Hopefuly they were anyway!
    //
    k3 = 0;
    c2 = s2 = c3 = s3 = 0.0f;
    
	a = &(data[0].r)-1;
	b = &(data[0].i)-1;
    
    
    // the following two constants should agree with the array dimensions.
    //
	maxf=23;
	maxp=209;
	if(n < 2) goto l998;
	inc=isn;
	rad=8.0*atan(1.0);
	s72=rad/5.0;
	c72=cos(s72);
	s72=sin(s72);
	s120=sqrt(0.75);
	if(isn >= 0) goto l10;
	s72=-s72;
	s120=-s120;
	rad=-rad;
	inc=-inc;
l10:	nt=inc*ntot;
	ks=inc*nspan;
	kspan=ks;
	nn=nt-inc;
	jc=ks/n;
	radf=rad*jc*0.5;
	i=0;
	jf=0;
    //	determine the factors of n
    //
	m=0;
	k=n;
	goto l20;
l15:	m=m+1;
	nfac[m]=4;
	k=k/16;
l20:	if(k-(k/16)*16 == 0) goto l15;
	j=3;
	jj=9;
	goto l30;
l25:	m=m+1;
	nfac[m]=j;
	k=k/jj;
l30:	if(k%jj == 0) goto l25;
	j=j+2;
	jj=j*j;	// changed from jj=j**2 by a.blake 26/1/93
	if(jj <= k) goto l30;
	if(k > 4) goto l40;
	kt=m;
	nfac[m+1]=k;
	if(k != 1) m=m+1;
	goto l80;
l40:	if(k-(k/4)*4 != 0) goto l50;
	m=m+1;
	nfac[m]=2;
	k=k/4;
l50:	kt=m;
	j=2;
l60:	if(k%j != 0) goto l70;
	m=m+1;
	nfac[m]=j;
	k=k/j;
l70:	j=((j+1)/2)*2+1;
	if(j <= k) goto l60;
l80:	if(kt == 0) goto l100;
	j=kt;
l90:	m=m+1;
	nfac[m]=nfac[j];
	j=j-1;
	if(j != 0) goto l90;
    //	compute fourier transform
    //
l100:	sd=radf/( double )kspan;
	cd=sin(sd);
	cd=cd*cd*2.0;
	sd=sin(sd+sd);
	kk=1;
	i=i+1;
	if(nfac[i] != 2) goto l400;
    //	transform for factor of 2 (including rotation factor)
    //
	kspan=kspan/2;
	k1=kspan+2;
l210:	k2=kk+kspan;
	ak=a[k2];
	bk=b[k2];
	a[k2]=a[kk]-ak;
	b[k2]=b[kk]-bk;
	a[kk]=a[kk]+ak;
	b[kk]=b[kk]+bk;
	kk=k2+kspan;
	if(kk <= nn) goto l210;
	kk=kk-nn;
	if(kk <= jc) goto l210;
	if(kk > kspan) goto l800;
l220:	c1=1.0-cd;
	s1=sd;
l230:	k2=kk+kspan;
	ak=a[kk]-a[k2];
	bk=b[kk]-b[k2];
	a[kk]=a[kk]+a[k2];
	b[kk]=b[kk]+b[k2];
	a[k2]=c1*ak-s1*bk;
	b[k2]=s1*ak+c1*bk;
	kk=k2+kspan;
	if(kk < nt) goto l230;
	k2=kk-nt;
	c1=-c1;
	kk=k1-k2;
	if(kk > k2) goto l230;
	ak=c1-(cd*c1+sd*s1);
	s1=(sd*c1-cd*s1)+s1;
    
    // the following three statements compensate for truncation
    // error. if rounded arithmetic is used, substitute
    // c1=ak
    //
#ifdef ROUND
	c1=ak;
#else
	c1=0.5/(ak*ak+s1*s1)+0.5;
	s1=c1*s1;
	c1=c1*ak;
#endif
	kk=kk+jc;
	if(kk < k2) goto l230;
	k1=k1+inc+inc;
	kk=(k1-kspan)/2+jc;
	if(kk <= jc+jc) goto l220;
	goto l100;
    // transform for factor of 3 (optional code)
    //
l320:	k1=kk+kspan;
	k2=k1+kspan;
	ak=a[kk];
	bk=b[kk];
	aj=a[k1]+a[k2];
	bj=b[k1]+b[k2];
	a[kk]=ak+aj;
	b[kk]=bk+bj;
	ak=-0.5*aj+ak;
	bk=-0.5*bj+bk;
	aj=(a[k1]-a[k2])*s120;
	bj=(b[k1]-b[k2])*s120;
	a[k1]=ak-bj;
	b[k1]=bk+aj;
	a[k2]=ak+bj;
	b[k2]=bk-aj;
	kk=k2+kspan;
	if(kk < nn) goto l320;
	kk=kk-nn;
	if(kk <= kspan) goto l320;
	goto l700;
    // transform for factor of 4
    //
l400:	if(nfac[i] != 4) goto l600;
	kspnn=kspan;
	kspan=kspan/4;
l410:	c1=1.0;
	s1=0.0;
l420:	k1=kk+kspan;
	k2=k1+kspan;
	k3=k2+kspan;
	akp=a[kk]+a[k2];
	akm=a[kk]-a[k2];
	ajp=a[k1]+a[k3];
	ajm=a[k1]-a[k3];
	a[kk]=akp+ajp;
	ajp=akp-ajp;
	bkp=b[kk]+b[k2];
	bkm=b[kk]-b[k2];
	bjp=b[k1]+b[k3];
	bjm=b[k1]-b[k3];
	b[kk]=bkp+bjp;
	bjp=bkp-bjp;
	if(isn < 0) goto l450;
	akp=akm-bjm;
	akm=akm+bjm;
	bkp=bkm+ajm;
	bkm=bkm-ajm;
	if(s1 == 0.) goto l460;
l430:	a[k1]=akp*c1-bkp*s1;
	b[k1]=akp*s1+bkp*c1;
	a[k2]=ajp*c2-bjp*s2;
	b[k2]=ajp*s2+bjp*c2;
	a[k3]=akm*c3-bkm*s3;
	b[k3]=akm*s3+bkm*c3;
	kk=k3+kspan;
	if(kk <= nt) goto l420;
l440:	c2=c1-(cd*c1+sd*s1);
	s1=(sd*c1-cd*s1)+s1;
    // the following three statements compensate for truncation
    // error. if rounded arithmetic is used, substitute
    // c1=c2
    //
#ifdef ROUND
	c1=c2;
#else
	c1=0.5/(c2*c2+s1*s1)+0.5;
	s1=c1*s1;
	c1=c1*c2;
#endif
	c2=c1*c1-s1*s1;
	s2=2.0*c1*s1;
	c3=c2*c1-s2*s1;
	s3=c2*s1+s2*c1;
	kk=kk-nt+jc;
	if(kk <= kspan) goto l420;
	kk=kk-kspan+inc;
	if(kk <= jc) goto l410;
	if(kspan == jc) goto l800;
	goto l100;
l450:	akp=akm+bjm;
	akm=akm-bjm;
	bkp=bkm-ajm;
	bkm=bkm+ajm;
	if(s1 != 0.) goto l430;
l460:	a[k1]=akp;
	b[k1]=bkp;
	a[k2]=ajp;
	b[k2]=bjp;
	a[k3]=akm;
	b[k3]=bkm;
	kk=k3+kspan;
	if(kk <= nt) goto l420;
	goto l440;
    // transform for factor of 5 (optional code)
    //
l510:	c2=c72*c72-s72*s72;
	s2=2.0*c72*s72;
l520:	k1=kk+kspan;
	k2=k1+kspan;
	k3=k2+kspan;
	k4=k3+kspan;
	akp=a[k1]+a[k4];
	akm=a[k1]-a[k4];
	bkp=b[k1]+b[k4];
	bkm=b[k1]-b[k4];
	bkp=b[k1]+b[k4];
	bkm=b[k1]-b[k4];
	ajp=a[k2]+a[k3];
	ajm=a[k2]-a[k3];
	bjp=b[k2]+b[k3];
	bjm=b[k2]-b[k3];
	aa=a[kk];
	bb=b[kk];
	a[kk]=aa+akp+ajp;
	b[kk]=bb+bkp+bjp;
	ak=akp*c72+ajp*c2+aa;
	bk=bkp*c72+bjp*c2+bb;
	aj=akm*s72+ajm*s2;
	bj=bkm*s72+bjm*s2;
	a[k1]=ak-bj;
	a[k4]=ak+bj;
	b[k1]=bk+aj;
	b[k4]=bk-aj;
	ak=akp*c2+ajp*c72+aa;
	bk=bkp*c2+bjp*c72+bb;
	aj=akm*s2-ajm*s72;
	bj=bkm*s2-bjm*s72;
	a[k2]=ak-bj;
	a[k3]=ak+bj;
	b[k2]=bk+aj;
	b[k3]=bk-aj;
	kk=k4+kspan;
	if(kk < nn) goto l520;
	kk=kk-nn;
	if(kk <= kspan) goto l520; else goto l700;
    // transform for odd factors
    //
l600:	k=nfac[i];
	kspnn=kspan;
	kspan=kspan/k;
	if(k == 3) goto l320;
	if(k == 5) goto l510;
	if(k == jf) goto l640;
	jf=k;
	s1=rad/( double ) k;
	c1=cos(s1);
	s1=sin(s1);
	if(jf > maxf) goto l998;
	ck[jf]=1.0;
	sk[jf]=0.;
	j=1;
l630:	ck[j]=ck[k]*c1+sk[k]*s1;
	sk[j]=ck[k]*s1-sk[k]*c1;
	k=k-1;
	ck[k]=ck[j];
	sk[k]=-sk[j];
	j=j+1;
	if(j < k) goto l630;
l640:	k1=kk;
	k2=kk+kspnn;
	aa=a[kk];
	bb=b[kk];
	ak=aa;
	bk=bb;
	j=1;
	k1=k1+kspan;
l650:	k2=k2-kspan;
	j=j+1;
	at[j]=a[k1]+a[k2];
	ak=at[j]+ak;
	bt[j]=b[k1]+b[k2];
	bk=bt[j]+bk;
	j=j+1;
	at[j]=a[k1]-a[k2];
	bt[j]=b[k1]-b[k2];
	k1=k1+kspan;
	if(k1 < k2) goto l650;
	a[kk]=ak;
	b[kk]=bk;
	k1=kk;
	k2=kk+kspnn;
	j=1;
l660:	k1=k1+kspan;
	k2=k2-kspan;
	jj=j;
	ak=aa;
	bk=bb;
	aj=0.;
	bj=0.;
	k=1;
l670:	k=k+1;
	ak=at[k]*ck[jj]+ak;
	bk=bt[k]*ck[jj]+bk;
	k=k+1;
	aj=at[k]*sk[jj]+aj;
	bj=bt[k]*sk[jj]+bj;
	jj=jj+j;
	if(jj > jf) jj=jj-jf;
	if(k < jf) goto l670;
	k=jf-j;
	a[k1]=ak-bj;
	b[k1]=bk+aj;
	a[k2]=ak+bj;
	b[k2]=bk-aj;
	j=j+1;
	if(j < k) goto l660;
	kk=kk+kspnn;
	if(kk <= nn) goto l640;
	kk=kk-nn;
	if(kk <= kspan) goto l640;
    // multiply by rotation factor (except for factors of 2 and 4)
    //
l700:	if(i == m) goto l800;
	kk=jc+1;
l710:	c2=1.0-cd;
	s1=sd;
l720:	c1=c2;
	s2=s1;
	kk=kk+kspan;
l730:	ak=a[kk];
	a[kk]=c2*ak-s2*b[kk];
	b[kk]=s2*ak+c2*b[kk];
	kk=kk+kspnn;
	if(kk <= nt) goto l730;
	ak=s1*s2;
	s2=s1*c2+c1*s2;
	c2=c1*c2-ak;
	kk=kk-nt+kspan;
	if(kk <= kspnn) goto l730;
	c2=c1-(cd*c1+sd*s1);
	s1=s1+(sd*c1-cd*s1);
    // the following three statements compensate for truncation
    // error. if rounded arithmetic is used, they may
    // be deleted.
    // 
#ifndef ROUND
	c1=0.5/(c2*c2+s1*s1)+0.5;
	s1=c1*s1;
	c2=c1*c2;
#endif
	kk=kk-kspnn+jc;
	if(kk <= kspan) goto l720;
	kk=kk-kspan+jc+inc;
	if(kk <= jc+jc) goto l710;
	goto l100;
    // permute the results to normal order---done in two stages
    // permutation for square factors of n
    //
l800:	np[1]=ks;
	if(kt == 0) goto l890;
	k=kt+kt+1;
	if(m < k) k=k-1;
	j=1;
	np[k+1]=jc;
l810:	np[j+1]=np[j]/nfac[j];
	np[k]=np[k+1]*nfac[j];
	j=j+1;
	k=k-1;
	if(j < k) goto l810;
	k3=np[k+1];
	kspan=np[2];
	kk=jc+1;
	j=1;
	k2=kspan+1;
	if(n != ntot) goto l850;
    // permutation for single-variate transform (optional code)
    //
l820:	ak=a[kk];
	a[kk]=a[k2];
	a[k2]=ak;
	bk=b[kk];
	b[kk]=b[k2];
	b[k2]=bk;
	kk=kk+inc;
	k2=kspan+k2;
	if(k2 < ks) goto l820;
l830:	k2=k2-np[j];
	j=j+1;
	k2=np[j+1]+k2;
	if(k2 > np[j]) goto l830;
	j=1;
l840:	if(kk < k2) goto l820;
	kk=kk+inc;
	k2=kspan+k2;
	if(k2 < ks) goto l840;
	if(kk < ks) goto l830;
	jc=k3;
	goto l890;
    // permutation for multivariate transform
    //
l850:	k=kk+jc;
l860:	ak=a[kk];
	a[kk]=a[k2];
	a[k2]=ak;
	bk=b[kk];
	b[kk]=b[k2];
	b[k2]=bk;
	kk=kk+inc;
	k2=k2+inc;
	if(kk < k) goto l860;
	kk=kk+ks-jc;
	k2=k2+ks-jc;
	if(kk < nt) goto l850;
	k2=k2-nt+kspan;
	kk=kk-nt+jc;
	if(k2 < ks) goto l850;
l870:	k2=k2-np[j];
	j=j+1;
	k2=np[j+1]+k2;
	if(k2 > np[j]) goto l870;
	j=1;
l880:	if(kk < k2) goto l850;
	kk=kk+jc;
	k2=kspan+k2;
	if(k2 < ks) goto l880;
	if(kk < ks) goto l870;
	jc=k3;
l890:	if(2*kt+1 >= m) return NO_ERROR;
	kspnn=np[kt+1];
    // permutation for square-free factors of n
    //
	j=m-kt;
	nfac[j+1]=1;
l900:	nfac[j]=nfac[j]*nfac[j+1];
	j=j-1;
	if(j != kt) goto l900;
	kt=kt+1;
	nn=nfac[kt]-1;
	if(nn > maxp) goto l998;
	jj=0;
	j=0;
	goto l906;
l902:	jj=jj-k2;
	k2=kk;
	k=k+1;
	kk=nfac[k];
l904:	jj=kk+jj;
	if(jj >= k2) goto l902;
	np[j]=jj;
l906:	k2=nfac[kt];
	k=kt+1;
	kk=nfac[k];
	j=j+1;
	if(j <= nn) goto l904;
    // determine the permutation cycles of length greater than 1
    //
	j=0;
	goto l914;
l910:	k=kk;
	kk=np[k];
	np[k]=-kk;
	if(kk != j) goto l910;
	k3=kk;
l914:	j=j+1;
	kk=np[j];
	if(kk < 0) goto l914;
	if(kk != j) goto l910;
	np[j]=-j;
	if(j != nn) goto l914;
	maxf=inc*maxf;
    // reorder a and b, following the permutation cycles
    //
	goto l950;
l924:	j=j-1;
	if(np[j] < 0) goto l924;
	jj=jc;
l926:	kspan=jj;
	if(jj > maxf) kspan=maxf;
	jj=jj-kspan;
	k=np[j];
	kk=jc*k+i+jj;
	k1=kk+kspan;
	k2=0;
l928:	k2=k2+1;
	at[k2]=a[k1];
	bt[k2]=b[k1];
	k1=k1-inc;
	if(k1 != kk) goto l928;
l932:	k1=kk+kspan;
	k2=k1-jc*(k+np[k]);
	k=-np[k];
l936:	a[k1]=a[k2];
	b[k1]=b[k2];
	k1=k1-inc;
	k2=k2-inc;
	if(k1 != kk) goto l936;
	kk=k2;
	if(k != j) goto l932;
	k1=kk+kspan;
	k2=0;
l940:	k2=k2+1;
	a[k1]=at[k2];
	b[k1]=bt[k2];
	k1=k1-inc;
	if(k1 != kk) goto l940;
	if(jj != 0) goto l926;
	if(j != 1) goto l924;
l950:	j=k3+1;
	nt=nt-kspnn;
	i=nt-inc+1;
	if(nt >= 0) goto l924;
	return NO_ERROR;
    // error finish, insufficient array storage
    //
l998:	isn=0;
	return INVALID_FFT_SIZE;
}

