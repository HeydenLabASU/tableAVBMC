#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double vec[3];

typedef int idx[3];

typedef struct {
    int restart;
	int cycles;
	double step;
	double bias;
	double r1,r1Sq;
	double r2,r2Sq;
	double rex,rexSq;
	double temp;
	double Vin,Vout;
	double ratioAVBMC;
	int printFreq;
	char fnXYZ[300];
} t_mc;

typedef struct {
        double rmin;
        double delta;
        int n;
        double *pot;
} t_pot;

typedef struct {
	int nPtot;
	int nTypes;
	int nP[10];
	double conc;
	double box;
	vec *crd;
	int *type;
	t_pot **pot;
	double rcut,rcutSq;
} t_sys;

typedef struct {
	int *list;
	int n;
} t_list;

typedef struct {
	int nGrid,nGrid3;
	double dGrid;
	idx *gridIdx;
	t_list *gridList;
} t_grid;

int fatal(const char *info)
{
        FILE *out;
        out=fopen("FATAL_ERROR.LOG","w");
        fprintf(out,"%s\n",info);
        fclose(out);
        printf("%s\n",info);
        exit(1);
}

int saveOpenRead(FILE **io,char *fn)
{
        char buffer[200];
        if((io[0]=fopen(fn,"r"))==NULL)
        {
                sprintf(buffer,"FILE NOT FOUND!!!: %s\n -> saveOpenRead\n -> io.c\n",fn);
                fatal(buffer);
        }
        return 0;
}

void *save_malloc(size_t size)
{
        void *ptr;
        char buffer[200];

        if((ptr=malloc(size))==NULL)
        {
                sprintf(buffer,"MEMORY ALLOCATION FAILED!!!\ntried to allocate %lu bytes\n -> save_malloc\n -> alloc.c\n",size);
                fatal(buffer);
        }

        return ptr;
}

int getLineFromCOM(FILE *in,char *buffer,int max)
{
        if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        while(strncmp(buffer,"#",1)==0)
        {
                if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        }
        return 0;
}

int readPot(char *fn,t_pot *pot) {
	FILE *io;
	char buffer[300];
	float tmp1,tmp2;
	int i;

	saveOpenRead(&io,fn);
	pot[0].n=0;
	while(fgets(buffer,300,io)!=NULL) {
		pot[0].n++;
	}
	rewind(io);
	pot[0].pot=(double*)save_malloc(pot[0].n*sizeof(double));
	for(i=0;i<pot[0].n;i++) {
		fgets(buffer,300,io);
		sscanf(buffer,"%f %f",&tmp1,&tmp2);
		if(i==0) pot[0].rmin=(double)tmp1;
		if(i==1) pot[0].delta=((double)tmp1)-pot[0].rmin;
		pot[0].pot[i]=(double)tmp2;
	}
	fclose(io);
}

int readXYZ(char *fn,t_sys *sys) {
	FILE *io;
	char buffer[300];
	int i,j;
	char pName[10];
	char pNameCur[10];
	float x,y,z;

	saveOpenRead(&io,fn);
	fgets(buffer,300,io);
	sscanf(buffer,"%d",&sys[0].nPtot);
	sys[0].crd=(vec*)save_malloc(sys[0].nPtot*sizeof(vec));
	fgets(buffer,300,io);

	sys[0].nTypes=0;
	j=-1;
	for(i=0;i<sys[0].nPtot;i++) {
		fgets(buffer,300,io);
		sscanf(buffer,"%s %f %f %f",pName,&x,&y,&z);
		sys[0].crd[i][0]=(double)x;
		sys[0].crd[i][1]=(double)y;
		sys[0].crd[i][2]=(double)z;
		if(i==0 || strcmp(pName,pNameCur)!=0) {
			if(sys[0].nTypes==10) fatal("only <=10 distinct particle types supported\n");
			sys[0].nTypes++;
			j++;
			sys[0].nP[j]=1;
			strcpy(pNameCur,pName);
		} else {
			sys[0].nP[j]++;
		}
	}
	fclose(io);

	printf("read XYZ file: %s\n",fn);
	for(i=0;i<sys[0].nTypes;i++) {
		printf(" - found %d particles of type %d\n",sys[0].nP[i],i+1);
	}
	printf("\n");
	return 0;
}
	
int getInput(char *fnCOM,t_sys *sys,t_mc *mc,double *maxMem) {
        FILE *io;
	int restart;
        char buffer[300];
        float tmp;
	int i=1;
	int j,k;
	char fnXYZ[300];
	char fnPot[300];
	int n,x,y,z;
	double d;

        saveOpenRead(&io,fnCOM);

	getLineFromCOM(io,buffer,300);
	sscanf(buffer,"%d",&restart);
	printf("%2d -> read %20s : %d\n",i,"restart",restart);i++;
    mc[0].restart=restart;

	getLineFromCOM(io,buffer,300);
	sscanf(buffer,"%s",fnXYZ);
	printf("%2d -> read %20s : %s\n",i,"fnXYZ",fnXYZ);i++;
	if(restart==1) readXYZ(fnXYZ,sys);

	getLineFromCOM(io,buffer,300);
	if(restart==0) {
		sscanf(buffer,"%d",&sys[0].nTypes);
		printf("%2d -> read %20s : %d\n",i,"nTypes",sys[0].nTypes);i++;
		if(sys[0].nTypes>10) fatal("only <=10 distinct particle types supported\n");
		sys[0].nPtot=0;
	}
	for(j=0;j<sys[0].nTypes;j++) {
		getLineFromCOM(io,buffer,300);
		if(restart==0) {
			sscanf(buffer,"%d",&sys[0].nP[j]);
			printf("%2d -> read %17s[%d] : %d\n",i,"nP",j,sys[0].nP[j]);i++;
			sys[0].nPtot+=sys[0].nP[j];
		}
	}

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
	printf("%2d -> read %20s : %f\n",i,"conc [mol/L]",tmp);i++;
        sys[0].conc=(double)tmp;
	sys[0].box=pow(sys[0].nPtot/(sys[0].conc*6.022*0.0001),1.0/3.0);

	if(restart==0) {
		sys[0].crd=(vec*)save_malloc(sys[0].nPtot*sizeof(vec));
		n=1;
		while(n*n*n<sys[0].nPtot) n++;
		d=sys[0].box/n;
		j=0;
		for(x=0;x<n;x++) {
			for(y=0;y<n;y++) {
				for(z=0;z<n;z++) {
					if(j<sys[0].nPtot) {
						sys[0].crd[j][0]=0.5*d+x*d;
						sys[0].crd[j][1]=0.5*d+y*d;
						sys[0].crd[j][2]=0.5*d+z*d;
					}
					j++;
				}
			}
		}
		/*random position of particles instead of regular grid*/
		/*for(j=0;j<sys[0].nPtot;j++) {
			sys[0].crd[j][0]=(((double)rand())/RAND_MAX)*sys[0].box;
			sys[0].crd[j][1]=(((double)rand())/RAND_MAX)*sys[0].box;
			sys[0].crd[j][2]=(((double)rand())/RAND_MAX)*sys[0].box;
		}*/
	}

	getLineFromCOM(io,buffer,300);
	sscanf(buffer,"%f",&tmp);
	printf("%2d -> read %20s : %f\n",i,"temperature [K]",tmp);i++;
	mc[0].temp=(double)tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"bias",tmp);i++;
        mc[0].bias=(double)tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"r1 [A]",tmp);i++;
        mc[0].r1=(double)tmp;
	mc[0].r1Sq=mc[0].r1*mc[0].r1;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"r2 [A]",tmp);i++;
        mc[0].r2=(double)tmp;
	mc[0].r2Sq=mc[0].r2*mc[0].r2;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"r_excl [A]",tmp);i++;
        mc[0].rex=(double)tmp;
	mc[0].rexSq=mc[0].rex*mc[0].rex;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"r_cutoff [A]",tmp);i++;
        sys[0].rcut=(double)tmp;
        sys[0].rcutSq=sys[0].rcut*sys[0].rcut;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"MC step size [A]",tmp);i++;
        mc[0].step=(double)tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",&mc[0].cycles);
        printf("%2d -> read %20s : %d\n",i,"MC cycles",mc[0].cycles);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"ratio of AVBMC moves",tmp);i++;
        mc[0].ratioAVBMC=(double)tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",&mc[0].printFreq);
        printf("%2d -> read %20s : %d\n",i,"print frequency",mc[0].printFreq);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",mc[0].fnXYZ);
        printf("%2d -> read %20s : %s\n",i,"xyz output file",mc[0].fnXYZ);i++;

	sys[0].pot=(t_pot**)save_malloc(sys[0].nTypes*sizeof(t_pot*));
	for(j=0;j<sys[0].nTypes;j++) {
		sys[0].pot[j]=(t_pot*)save_malloc(sys[0].nTypes*sizeof(t_pot));
	}
	for(j=0;j<sys[0].nTypes;j++) {
		for(k=j;k<sys[0].nTypes;k++) {
			getLineFromCOM(io,buffer,300);
			sscanf(buffer,"%s",fnPot);
			printf("%2d -> read %14s[%d][%d] : %s\n",i,"potential",j,k,fnPot);i++;
			readPot(fnPot,&sys[0].pot[j][k]);
			if(j!=k) readPot(fnPot,&sys[0].pot[k][j]);
		}
	}

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %e\n",i,"max memory [bytes]",tmp);i++;
        maxMem[0]=(double)tmp;
	fclose(io);

	sys[0].type=(int*)save_malloc(sys[0].nPtot*sizeof(int));
	k=0;
	for(i=0;i<sys[0].nTypes;i++) {
		for(j=0;j<sys[0].nP[i];j++) {
			sys[0].type[k]=i;
			k++;
		}
	}

	mc[0].Vin=4.0/3.0*3.14159*mc[0].r2*mc[0].r2*mc[0].r2;
	mc[0].Vin-=4.0/3.0*3.14159*mc[0].r1*mc[0].r1*mc[0].r1;
	mc[0].Vout=sys[0].box*sys[0].box*sys[0].box-mc[0].Vin;
	return 0;
}

int vecSub(vec a,vec b,vec *c) {
	c[0][0]=a[0]-b[0];
	c[0][1]=a[1]-b[1];
	c[0][2]=a[2]-b[2];
	return 0;
}

int vecAdd(vec a,vec b,vec *c) {
        c[0][0]=a[0]+b[0];
        c[0][1]=a[1]+b[1];
        c[0][2]=a[2]+b[2];
        return 0;
}

int vecCpy(vec a,vec *b) {
	b[0][0]=a[0];
	b[0][1]=a[1];
	b[0][2]=a[2];
	return 0;
}

int vecNorm(vec a,double *r) {
	r[0]=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	return 0;
}

int vecNormSq(vec a, double *rSq) {
	rSq[0]=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
	return 0;
}

int distPBC(vec a,vec b,double *r,double box) {
	vec link;
	int i;

	vecSub(a,b,&link);
	for(i=0;i<3;i++) {
                link[i]=link[i]-rint(link[i]/box)*box;
        }
        vecNorm(link,r);
	return 0;
}

int distSqPBC(vec a,vec b,double *rSq,double box) {
        vec link;
        int i;

        vecSub(a,b,&link);
        for(i=0;i<3;i++) {
                link[i]=link[i]-rint(link[i]/box)*box;
        }
        vecNormSq(link,rSq);
        return 0;
}

int eij(vec icrd,t_sys sys,t_mc mc,int i,int j,double *energy) {
	vec link;
	double r;
	int t1,t2;
	double rmin;
	double delta;
	int n;
	int k;
	double lambda;

	t1=sys.type[i];
	t2=sys.type[j];

	rmin=sys.pot[t1][t2].rmin;
	delta=sys.pot[t1][t2].delta;
	n=sys.pot[t1][t2].n;

	distPBC(icrd,sys.crd[j],&r,sys.box);
	if(r<mc.rex) {
		energy[0]=99999.99;
		return 1;
	} else {
		if(r>sys.rcut) energy[0]=0.0;
		else if(r<rmin) {
			energy[0]=sys.pot[t1][t2].pot[0];
		} else {
			k=(int)((r-rmin)/delta);
			if(k<n-1) {
				lambda=(r-(rmin+k*delta))/delta;
				energy[0]=lambda*sys.pot[t1][t2].pot[k];
				energy[0]+=(1-lambda)*sys.pot[t1][t2].pot[k+1];
			} else energy[0]=0.0;
		}
		return 0;
	}
}

int getGridIdx(vec a,double dGrid,idx *idx) {
	idx[0][0]=(int)(a[0]/dGrid);
	idx[0][1]=(int)(a[1]/dGrid);
	idx[0][2]=(int)(a[2]/dGrid);
	return 0;
}

int g2l(idx idx,int n) {
	return n*(idx[0]*n+idx[1])+idx[2];
}

int rmFromList(t_list *l,int i) {
	t_list tmp;
	int j,k;
	
	tmp.list=save_malloc(l[0].n*sizeof(int));
	tmp.n=l[0].n;
	
	for(j=0;j<l[0].n;j++) {
		tmp.list[j]=l[0].list[j];
	}
	k=0;
	for(j=0;j<tmp.n;j++) {
		if(tmp.list[j]!=i) {
			l[0].list[k]=tmp.list[j];
			k++;
		}
	}
	l[0].n=k;
	free(tmp.list);
	return 0;
}

int appendToList(t_list *l,int i) {
	int n;

	n=l[0].n;
	l[0].list[n]=i;
	l[0].n++;
	return 0;
}

int closeVoxel(idx center,idx *all,int nGrid) {
	int i,x,y,z;
	i=0;
        for(x=-1;x<=1;x++) {
                for(y=-1;y<=1;y++) {
                        for(z=-1;z<=1;z++) {
                                all[i][0]=center[0]+x;
                                if(all[i][0]<0) all[i][0]+=nGrid;
                                else if(all[i][0]>=nGrid) all[i][0]-=nGrid;

                                all[i][1]=center[1]+y;
                                if(all[i][1]<0) all[i][1]+=nGrid;
                                else if(all[i][1]>=nGrid) all[i][1]-=nGrid;

                                all[i][2]=center[2]+z;
                                if(all[i][2]<0) all[i][2]+=nGrid;
                                else if(all[i][2]>=nGrid) all[i][2]-=nGrid;

                                i++;
                        }
                }
        }
	return 0;
}

int getParticleEnergy(vec ipcrd,t_sys sys,t_mc mc,int ip,t_grid g,double *eSum) {
        idx center;
        idx all[27];
        double ener;
        int i,j,k,l;
	int clash;
	int anyClash=0;

	getGridIdx(ipcrd,g.dGrid,&center);
	closeVoxel(center,all,g.nGrid);

	eSum[0]=0.0;
        for(i=0;i<27;i++) {
                j=g2l(all[i],g.nGrid);
                for(k=0;k<g.gridList[j].n;k++) {
                        l=g.gridList[j].list[k];
			if(ip!=l) {
                        	clash=eij(ipcrd,sys,mc,ip,l,&ener);
				if(clash==1) {
					anyClash++;
				}
				eSum[0]+=ener;
			}
                }
        }
        return anyClash;
}

int listInOut(t_sys sys,int ip,t_grid g,t_mc mc,t_list *in,t_list *out,int *tag) {
	idx center;
        idx all[27];
	int i,k,j,l;
	double rSq;

	for(i=0;i<sys.nPtot;i++) tag[i]=0;
	getGridIdx(sys.crd[ip],g.dGrid,&center);
	closeVoxel(center,all,g.nGrid);

	in[0].n=0;
	out[0].n=0;
	for(i=0;i<27;i++) {
                j=g2l(all[i],g.nGrid);
		for(k=0;k<g.gridList[j].n;k++) {
			l=g.gridList[j].list[k];
			if(ip!=l) {
				distSqPBC(sys.crd[ip],sys.crd[l],&rSq,sys.box);
				if(rSq>=mc.r1Sq && rSq<mc.r2Sq) {
					in[0].list[in[0].n]=l;
					in[0].n++;
					tag[l]=1;
				}
			}
		}
	}
	for(i=0;i<sys.nPtot;i++) {
		if(tag[i]==0) {
			out[0].list[out[0].n]=i;
			out[0].n++;
		}
	}
	return 0;
}

int randVec(double l,vec *vec) {
	double len=2.0;
	while(len>1.0) {
		vec[0][0]=2*(((double)rand())/RAND_MAX)-1.0;
		vec[0][1]=2*(((double)rand())/RAND_MAX)-1.0;
		vec[0][2]=2*(((double)rand())/RAND_MAX)-1.0;
		vecNormSq(vec[0],&len);
	}
	vec[0][0]*=l;
	vec[0][1]*=l;
	vec[0][2]*=l;
	return 0;
}

int placeInBox(vec *a,double box) {
	int i;
	for(i=0;i<3;i++) {
		if(a[0][i]<0) a[0][i]+=box;
		else if(a[0][i]>=box) a[0][i]-=box;
	}
	return 0;
}

int randIntRange(int max) {
	int tmp;
	tmp=((double)rand()/(RAND_MAX))*max;
	if(tmp==max) tmp--;
	return tmp;
}

int update(vec trial,int ip,t_sys sys,t_grid g) {
	int iOld,iNew;

	vecCpy(trial,&sys.crd[ip]);

	iOld=g2l(g.gridIdx[ip],g.nGrid);
	getGridIdx(trial,g.dGrid,&g.gridIdx[ip]);
	iNew=g2l(g.gridIdx[ip],g.nGrid);

	if(iOld!=iNew) {
		rmFromList(&g.gridList[iOld],ip);
		appendToList(&g.gridList[iNew],ip);
	}
	return 0;
}

int MCmove(t_sys sys,t_mc mc,t_grid g) {
	vec disp,trial;
	int sel;
	double eOld,eNew;
	double prob,dice;
	int i;
	int accepted=0;
	int clashOld,clashNew;
	
	sel=randIntRange(sys.nPtot);

	randVec(mc.step,&disp);
	vecAdd(sys.crd[sel],disp,&trial);
	placeInBox(&trial,sys.box);

	clashOld=getParticleEnergy(sys.crd[sel],sys,mc,sel,g,&eOld);
	clashNew=getParticleEnergy(trial,sys,mc,sel,g,&eNew);

	if(clashOld>clashNew) prob=1.0;
	else if(clashOld<clashNew) prob=0.0;
	else if(eNew<=eOld) prob=1.0;
	else prob=exp((eOld-eNew)/(0.0083145*mc.temp));

	dice=((double)rand())/RAND_MAX;
	if(dice<prob) {
		update(trial,sel,sys,g);
		accepted=1;
	}
	return accepted;
}

int AVBMCmove(t_sys sys,t_mc mc,t_grid g,t_list *in,t_list *out,int *tag) {
	int target,sel,tmp;
	double dice;
	double dSq;
	vec disp;
	vec trial;
	double eOld,eNew;
	double prob;
	int isOut;
	int accepted=0;
	int clashOld,clashNew;

	target=randIntRange(sys.nPtot);

	listInOut(sys,target,g,mc,in,out,tag);
	dice=((double)rand())/RAND_MAX;
	if(dice<mc.bias) { /*out->in*/
		if(out[0].n!=0) {
			tmp=randIntRange(out[0].n);
			sel=out[0].list[tmp];
			dSq=-1.0;
			while(dSq<=mc.r1Sq) {
				randVec(mc.r2,&disp);
				vecNormSq(disp,&dSq);
			}
			vecAdd(sys.crd[target],disp,&trial);
			placeInBox(&trial,sys.box);
			clashOld=getParticleEnergy(sys.crd[sel],sys,mc,sel,g,&eOld);
			clashNew=getParticleEnergy(trial,sys,mc,sel,g,&eNew);
			if(clashOld>clashNew) prob=1.0;
			else if(clashOld<clashNew) prob=0.0;
			else {
				prob=((1-mc.bias)*mc.Vin*out[0].n*exp((eOld-eNew)/(0.0083145*mc.temp)));
				prob/=(mc.bias*mc.Vout*(in[0].n+1));
			}
		} else prob=0.0;
	} else { /*in->out*/
		if(in[0].n!=0.0) {
			tmp=randIntRange(in[0].n);
			sel=in[0].list[tmp];
			isOut=0;
			while(isOut==0) {
				trial[0]=(((double)rand())/RAND_MAX)*sys.box;
				trial[1]=(((double)rand())/RAND_MAX)*sys.box;
				trial[2]=(((double)rand())/RAND_MAX)*sys.box;
				distSqPBC(sys.crd[target],trial,&dSq,sys.box);
				if(dSq>mc.r2Sq || dSq<mc.r1Sq) isOut=1;
			}
			placeInBox(&trial,sys.box);
			clashOld=getParticleEnergy(sys.crd[sel],sys,mc,sel,g,&eOld);
			clashNew=getParticleEnergy(trial,sys,mc,sel,g,&eNew);
			if(clashOld>clashNew) prob=1.0;
			else if(clashOld<clashNew) prob=0.0;
			else {
				prob=(mc.bias*mc.Vout*in[0].n*exp((eOld-eNew)/(0.0083145*mc.temp)));
				prob/=((1-mc.bias)*mc.Vin*(out[0].n+1));
			}
		} else prob=0.0;
	}
	dice=((double)rand())/RAND_MAX;
	if(dice<=prob) {
		update(trial,sel,sys,g);
		accepted=1;
	}
	return accepted;
}

int writeXYZ(FILE *outXYZ,t_sys sys,int i,char *types) {
	int j,k,l;

	fprintf(outXYZ,"%d\nMC cycle %5d box %15.8e\n",sys.nPtot,i,sys.box);
	l=0;
	for(j=0;j<sys.nTypes;j++) {
		for(k=0;k<sys.nP[j];k++) {
			fprintf(outXYZ,"%c %15.8e %15.8e %15.8e\n",
				types[j],sys.crd[l][0],sys.crd[l][1],sys.crd[l][2]);
			l++;
		}
	}
	fflush(outXYZ);
	return 0;
}

int main(int argc, char *argv[]) {
	t_sys sys;
	t_mc mc;
	t_grid g;
	int i,j,k,l,m,n,o,p;
	double eOld;
	t_list in;
	t_list out;
	int *tag;
	double maxMem;
	int MCtry,MCacc;
	int AVBMCtry,AVBMCacc;
	double dice;
	FILE *outXYZ;
	char types[]="ABCDEFGHIJ";
	char buffer[305];
	
	if(argc!=2) fatal("provide simulation parameter input file\n");

	getInput(argv[1],&sys,&mc,&maxMem);

	in.list=(int*)save_malloc(sys.nPtot*sizeof(int));
	out.list=(int*)save_malloc(sys.nPtot*sizeof(int));
	tag=(int*)save_malloc(sys.nPtot*sizeof(int));

	g.nGrid=(int)floor(sys.box/sys.rcut);
	g.nGrid3=g.nGrid*g.nGrid*g.nGrid;
	printf("initial guess for nb list grid: %d %d %d\n",g.nGrid,g.nGrid,g.nGrid);
	while((double)g.nGrid3*(double)sys.nPtot*4.0>maxMem) {
		g.nGrid--;
		g.nGrid3=g.nGrid*g.nGrid*g.nGrid;
	}
	printf("will use %d %d %d grid for nb lists\n",g.nGrid,g.nGrid,g.nGrid);
	fflush(stdout);
	g.dGrid=sys.box/g.nGrid;
	g.gridIdx=(idx*)save_malloc(sys.nPtot*sizeof(idx));
	g.gridList=(t_list*)save_malloc(g.nGrid3*sizeof(t_list));
	for(i=0;i<g.nGrid3;i++) {
		g.gridList[i].list=(int*)save_malloc(sys.nPtot*sizeof(int));
		g.gridList[i].n=0;
	}
	for(i=0;i<sys.nPtot;i++) {
		getGridIdx(sys.crd[i],g.dGrid,&g.gridIdx[i]);
		j=g2l(g.gridIdx[i],g.nGrid);
		appendToList(&g.gridList[j],i);
	}

	outXYZ=fopen(mc.fnXYZ,"w");
	if(mc.restart==0) {
        writeXYZ(outXYZ,sys,0,types);
	    printf("wrote frame at cycle: %5d\n",0);fflush(stdout);
    }
	MCtry=0;
	MCacc=0;
	AVBMCtry=0;
	AVBMCacc=0;
	for(i=0;i<mc.cycles;i++) {
		for(j=0;j<sys.nPtot;j++) {
			dice=((double)rand())/RAND_MAX;
			if(dice<mc.ratioAVBMC) {
				AVBMCtry++;
				AVBMCacc+=AVBMCmove(sys,mc,g,&in,&out,tag);
			} else {
				MCtry++;
				MCacc+=MCmove(sys,mc,g);
			}
		}
		if((i+1)%mc.printFreq==0) {
			writeXYZ(outXYZ,sys,i+1,types);
			printf("wrote frame at cycle: %5d MC:%9d/%-9d AVBMC:%9d/%-9d\n",
				i+1,MCacc,MCtry,AVBMCacc,AVBMCtry);
			fflush(stdout);
		}
	}
	fclose(outXYZ);
	
	sprintf(buffer,"last_%s",mc.fnXYZ);
	outXYZ=fopen(buffer,"w");
	writeXYZ(outXYZ,sys,i,types);
	printf("wrote final frame to %s\n",buffer);
	fflush(stdout);
	fclose(outXYZ);

	return 0;
}
