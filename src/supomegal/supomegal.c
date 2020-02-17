/* This code is distributed as part of the source-code package SWIP
/* that accompanies Pasquet and Bodet (2017). The package can be downloaded 
/* from the Geophysics source-code archive at http://software.seg.org/2017/0008
/* Use of this code is subject to acceptance of the terms and conditions
/* that can be found at http://software.seg.org/disclaimer2.txt
/* Copyright (c) 2017 by the Society of Exploration Geophysicists.
/* Reference: Pasquet, S., Bodet, L. (2017) SWIP: an integrated workflow for
/* surface-wave dispersion inversion and profiling, Geophysics, 82(6), 1-15,
/* doi: 10.1190/geo2016-0625.1.

/*LB : p-omega stack based on Herrmann, 2002, Computer Programs in Seismology */
/*     http://www.eas.slu.edu/People/RBHerrmann/ComputerPrograms.html */

/* DISERSION IN FOURIER DOMAIN 	 Herrmann CPS 2002	*/
/*              Adnand Bitri 1999			*/
/*            Modif Ludovic Bodet 2004                  */
/*            Modif Sylvain Pasquet 2014                  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "taup.h"
#include <signal.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TPAD_FAC 4 /* //2 a determiner*/
#define XPAD_FAC 2 /* //2 a determiner*/
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#define PAD_FAC 1.1
/*********************** self documentation **********************/
char *sdoc[] = {
"UNDER CONSTRUCTION                                                                        ",
" SUPOMEGAL - p-omega stack						",
"                                                                       ",
"    supomegal <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
" dt=tr.dt (from header) 	time sampling interval (secs)           ",
" nx=ntr   (counted from data)	number of horizontal samples (traces)	",
" nray=250 			number of velocity for p-omega stack  	",
" coord=0			=1 use keys gx-sx to compute offsets     ",
" xmin=1.0			offset on first trace  (m)         	",
" dx=1.0				horizontal sampling interval (m)	",
" 	OR xmin=tr.offset(0) and dx=tr.offset(i+1) - tr.offset(i)",
" xsca=1.0 			scaling factor (distance=distance/xsca)",
" tsca=1. !!!problem if <1 !!!		time scaling factor (time=time/xsca)",
" vmin=100.0			Minimum velocity for p-omega transform (s/m)	",
" vmax=5000			Maximum velocity for p-omega transform (s/m)	",
" fmax=100 maximum frequency Hz (user defined)				",
" fmin=0.0 minimum       						",
" flip=0			=1 					",
" verbose=0	verbose = 1 echoes information			",
" 								",
" (!!! if the source is on the seismogram left side, dx & xmin must be negative)",
" 								",
" tmpdir= 	 if non-empty, use the value as a directory path",
"		 prefix for storing temporary files; else if the",
"	         the CWP_TMPDIR environment variable is set use	",
"	         its value for the path; else use tmpfile()	",
"	UNDER CONSTRUCTION        						",
"	Ref:  Herrmann, 2002, Computer Programs in Seismology ",
"             http://www.eas.slu.edu/People/RBHerrmann/ComputerPrograms.html ",
NULL};

/**************** end self doc ********************************/

static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *fds;


segy tr;

int
main(int argc, char **argv)
{
	int ix,it,ir;		/* loop counters */
	int ntr;		/* number of input traces */
	int nt;			/* number of time samples */
	int nx;			/* number of horizontal samples */
	int ntfft;		/* transform length			*/
	int nfby2p1;		/* nfft/2 + 1				*/
	int option,dct,flip;		/* flag for requested opeartion */
	float dt,tsca;               /* Time sample interval */
        float dx,dist[8192],xsca;    /* horizontal sample interval */
	float coord;
	/*modifié ici mais j'utilise encore ancienne version compilée SU37*/
	float xmin;		/* offset on first trace */
        float vmin;             /* Minimum velocity for p-omega transform */
        float vmax;             /* Maximum velocity for p-omega transform */
	int nray;		/* number of velocity for p-omega stack */
	float fmin,fmax;	/* Min/Max frequencies for p-omega transform */
	int iom,nflow,nfupr,nifreq; /* frequency limits & counters (computation)*/
	float dc;		/* velocity sampling interval */
	float d1;		/* output sample interval in Hz		*/
	float omega,omegap,omegar;
	float pray,vit;
	float fac,ca;
	float **in_traces;	/* array[nx][nt] of input traces */
	register float *q,**ouray;
	register complex *cfac, *cwork0, *cwork1, *dat;	/* complex transformed traces	*/
	int verbose;		/* flag for echoing information */
	char *tmpdir;		/* directory path for tmp files */
	cwp_Bool istmpdir=cwp_false;/* true for user-given path */
	complex cfac1;

        /* hook up getpar to handle the parameters */
        initargs(argc,argv);
        requestdoc(1);

	if (!getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		err("you can't write in %s (or it doesn't exist)", tmpdir);

        /* get info from first trace */
        if (!gettr(&tr))  err("can't get first trace");
        nt = tr.ns;
        dt = (float) tr.dt/1000000.0;

        /* Store traces in tmpfile while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;
		if (verbose) warn("putting temporary files in %s", directory);
	}
	ntr = 0;
	ix=0;
	if (!getparfloat("tsca",&tsca))		tsca = 1.0;
	if (!getparfloat("xsca",&xsca))		xsca = 1.0;
	if (!getparfloat("xmin",&xmin))		xmin=1;
	if (!getparfloat("dx",&dx))		    dx =1;
	if (!getparfloat("coord",&coord))	coord =0;
        do {
                ++ntr;
                efwrite(&tr, 1, HDRBYTES, headerfp);
                efwrite(tr.data, FSIZE, nt, tracefp);
		if (!getparfloat("dx",&dx))	{
			if (tr.offset == 0.0) { 
			dist[ix]=(xmin + 1.0*dx*ix)/xsca;
			}else 
			{	
			dist[ix]=fabs(tr.offset)/xsca;
		}
		}else {
			dist[ix]=(xmin + 1.0*dx*ix)/xsca;
		}
		if (ix > 0){
		if (!getparfloat("dx",&dx)) dx =abs(dist[ix] - dist[ix-1]);
		}
		if (coord==1) { 
		dist[ix]=sqrt(pow(tr.gx-tr.sx,2)+pow(tr.gy-tr.sy,2)+pow(tr.gelev-tr.selev,2))/xsca;
		}
		++ix;
        } while (gettr(&tr));
	if (!getparfloat("xmin",&xmin))		xmin=dist[0];

	/* get general flags and parameters and set defaults */
        if (!getparint("nx",&nx))          	nx = ntr;
        if (!getparint("option",&option))   option = 1;
        if (!getparfloat("vmin",&vmin))		vmin = 100.0;
        if (!getparfloat("vmax",&vmax))		vmax = 5000.0;
        if (!getparfloat("fmax",&fmax))		fmax = 100.0;
	if (!getparfloat("fmin",&fmin))		fmin = 0.0;
	if (!getparfloat("dt",&dt))			dt = dt;
        if (!getparint("nray",&nray))       nray = 250;
	if (!getparint("flip",&flip))       flip = 0;
	if (dt == 0.0)
		err("header field dt not set, must be getparred");
	dt = dt/tsca;

	/* Set up pfa fft */
	/*ntfft = npfa(nt);*/
	
	/*ntfft=npfao((int)(nt*PAD_FAC),2*nt);*/
	/*ntfft = npfaro(nt, LOOKFAC * nt);*/
	/*ntfft = npfar(nt);*/
	/* FFT */
	/*pfarc(1, nfft, rt, ct);*/

	if (ntfft >= SU_NFLTS || ntfft >= PFA_MAX)
		 err("Padded nt=%d--too big", ntfft);

	ntfft=pad2(nt);
	/*warn("%d%",ntfft);*/
	nfby2p1 = ntfft/2 + 1;
	d1 = 1.0/(ntfft*dt);
	nflow = (int)(fmin/d1+1);
	nfupr = (int)(fmax/d1+1+0.5);
	nifreq = nfupr-nflow +1;


	/* allocate space */
        in_traces = alloc2float(nt, ntr);
	cwork0 = ealloc1complex(ntfft);
	cwork1 = ealloc1complex(ntfft);
	cfac = ealloc1complex(ntfft/2);
	q = ealloc1float(nifreq+3);
	ouray = ealloc2float(nray,nifreq+3);
	dat = ealloc1complex(ntfft/2);

        /* load traces into an array and close temp file */
	erewind(headerfp);
        erewind(tracefp);
        for (ix=0; ix<ntr; ix++)
                fread (in_traces[ix], FSIZE, nt, tracefp);
        efclose (tracefp);
	if (istmpdir) eremove(tracefile);

	fds = fopen(".speclud","w");

	for (ix=0; ix<nx; ix++) {
		for(it=0;it<nt;it++) cwork1[it]=cmplx(in_traces[ix][it],0.0);
		/*pad t dimension with zeros*/
	        for(it=nt; it<ntfft ; it++) cwork1[it]=cmplx(0.0,0.0);
		/* FFT */
		/*pfacc(-1,ntfft,cwork1);*/
		frk(ntfft,cwork1,-1.0);
			fac =sqrt((double)(fabs(dist[ix])));
		for (it=0; it<ntfft; it++)cwork1[it]=cmul(cwork1[it],cmplx(fac,0.0));
		fwrite((char *) cwork1,1,ntfft*8,fds);
	}
	fclose(fds);

	fds = fopen(".speclud","r");
	/*warn("%f %f %f",dist[0], dist[1], dist[nx-1]);*/
	dc = (vmax - vmin)/(nray -1);
	dct=(int)(1000.0*dc);
	if (dct>=65535) warn("nray too small");
	for (ir=0; ir < nray; ir++) {
	rewind(fds);

		vit = vmin + ir*dc;
		pray = 1./vit;
		for (ix=0; ix<nx ; ix++){
			fread(cwork1,8,ntfft ,fds);
			for (iom=nflow;iom<nfupr; iom++){
				omega = 2.0*PI*iom/(ntfft*dt);
				omegap = omega*pray;	
				/*warn("%f %f %f",dist[0], dist[ix], dist[nx-1]);*/
				omegar = omegap*dist[ix];
				/*if (ix==0) cwork0[iom]=cwork1[iom];*/
				ca = rcabs(cwork1[iom]);
				if ( ca > 0.0)
					cfac1 = cdiv(cwork1[iom],cmplx(ca,0.0));
				else 
					cfac1 = cmplx(0.0,0.0);
				/*cfac1 = cmul(cfac1,cmplx((float)cos((double)omegar),(float)sin((double)omegar)));*/
				cfac1 = cmul(cfac1,cmplx(cos(omegar),sin(omegar)));
				if ( (float)(ix) == 0.0) {
	        			dat[iom] = cmplx(0.0,0.0);
					cfac[iom] = conjg(cfac1);
					}
				dat[iom] = cadd(dat[iom],cmul(cfac1,cfac[iom]));
				q[iom-nflow] = rcabs(dat[iom]);
				ouray[iom-nflow][ir] = rcabs(dat[iom]);
			}
		}
	if (flip==0) {
	     tr.dt=0;
	     tr.ns=nifreq;
	     tr.tracf=tr.tracr=tr.tracl=1+ir;
	     tr.d1 = d1*1.0;
	     tr.ntr = nray;
	     tr.f1 = fmin*1.0;
	     tr.d2 = dc*1.0;
	     tr.f2 = vmin*1.0;
		for (iom=0;iom<nifreq; iom++){
		tr.data[iom]=q[iom];
		}
		puttr(&tr);
	}
	}
	if (flip==1) {
	for (iom=0;iom<nifreq; iom++){
		/*tr.dt = (int)(dc*1000.0);*/
		tr.dt = dct;
		tr.ns=nray;
		tr.tracf=tr.tracr=tr.tracl=1+iom;
		/*tr.d1 = dc*1.0;*/
		tr.ntr = nifreq;
		/*tr.f1 =  vmin*1.0;*/
		tr.delrt =vmin*1.0;
		tr.d2 = d1*1.0;
		tr.f2 = fmin*1.0;
		for (ir=0; ir < nray; ir++) {
		tr.data[ir]=ouray[iom][ir];
		}
		puttr(&tr);
	}
	}
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);

	/* free allocated space */
	free2float(in_traces);
	system("rm -rf .speclud");

/*	return(CWP_Exit());*/

}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

frk(lx,cx,signi)
	complex cx[];
	float signi;
	int lx;
{
	complex car,cw,ct;	
	complex ctt;
	float rp,sc;
	float tt;
	int i,j,k,m;
	j=0;
	sc=sqrt(1./lx);
	for(i=0;i<lx;++i){
		if(i<=j){
			ct=cmul(cx[j],cmplx(sc,0.0));
			cx[j]=cmul(cx[i],cmplx(sc,0.0));
			cx[i]=ct;
			}
		m=lx/2;
		while(j>=m) {
			j-=m;
			m=m/2;
			if(m<1) break;
			}
		j+=m;
		}
	for(k=1;k<=(lx/2);k*=2){
		for(m=0;m<k;++m){
			rp=PI*signi*m/k;
			car=cmplx(0.0,rp);
			cw=cwp_cexp(car);
			for(i=m;i<lx;i+=(2*k)){
				ct=cmul(cw,cx[i+k]);
				cx[i+k]=csub(cx[i],ct);
				cx[i]=cadd(cx[i],ct);
				}
			}
		}
} 


int pad2(n)
	int n;
	{
	int pad;
	pad=1;
	while (pad< n){
		pad=pad*2;
		}
		return(pad);
	}




