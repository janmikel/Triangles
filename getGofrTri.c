/*
 CHANGE LOG
 
    May 21,2015 -Started this program
    
 OBJECTIVE:
    
    get g(r) for triangle (thiophene) system

 NOTES:
    1. Works only for 1 snapshot and a tape.
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

void alloc_r(float *& rr, int nAtom, int *& atype){
    int i;
    rr= (float *) malloc (nAtom*3 * sizeof(float));
    for(i=0;i<nAtom*3;i++) rr[i]=0;
    
    atype=(int *) malloc (nAtom * sizeof(int));
}
void alloc_bond(int *& bb, int nBond){
    int i;
    bb= (int *) malloc (nBond*2 * sizeof(int));
    for(i=0;i<nBond*2;i++) bb[i]=0;
}

int read_dump(char *dumpFile, float *& rr, float *& rru, float *pbcl, int *& atype){
    
    FILE *fpr0;
    int nAtom;
    int i,j,k;
    char buffer[256];

    fpr0=fopen(dumpFile,"r");
    
    while(1){
        fgets(buffer,254,fpr0);
        if(feof(fpr0)!=0){
            break;
        }
        else{
            for(i=0;i<3;i++) fgets(buffer,254,fpr0);
            sscanf(buffer,"%i",&nAtom);
            
            alloc_r(rr,nAtom,atype);   
            alloc_r(rru,nAtom,atype);
            
            fgets(buffer,254,fpr0);
            for(i=0;i<3;i++){ 
                fgets(buffer,254,fpr0);
                sscanf(buffer,"%f %f",&pbcl[i],&pbcl[i]);
            }
            fgets(buffer,254,fpr0);
            
            
            for(i=0;i<nAtom;i++){
                fgets(buffer,254,fpr0);
                int ID,at;
                float p[3],r[3];

                sscanf(buffer,"%i %i %f %f %f %f %f %f",&ID,&at,&p[0],&p[1],&p[2],&r[0],&r[1],&r[2]);
                //printf("Mark1\n");
                for(k=0;k<3;k++){ 
                    rr[(ID-1)*3+k]=p[k];
                    rru[(ID-1)*3+k]=r[k];
                }
                atype[(ID-1)]=at;
                
            }
            
        }
    }
    fclose(fpr0);
    return nAtom;
}

main(){
    int i,j,k;
    char buffer[256],TempC[30], AtomFile[256];
        
    while(1){
		fgets(buffer,254,stdin);
		if(feof(stdin)!=0){
			break;
		}
		else{
			if(strstr(buffer,"AtomFile")){
 				sscanf(buffer, "%s %s", &TempC,&AtomFile);
				printf("# %s %s\n",TempC,AtomFile);
			}
		}
    }
    
    int nAtom;
    float *rr0,*rru0;
    int *atype;
    float pbcl[3];
    int nMol;
    int NN = 6; //number of beads in one molecule
    
    nAtom=read_dump(AtomFile,rr0,rru0,pbcl,atype);
    printf("# pbcl: %f %f %f\n", pbcl[0],pbcl[1],pbcl[2]);
    printf("# nAtom: %i\n",nAtom);
    nMol = nAtom/NN;
    printf("# nMol: %i\n",nMol);
    
    float **rcm; //allocate center of mass of triangle
    rcm= (float**) malloc (nMol * sizeof(float*));
    for(i=0;i<nMol;i++){
        rcm[i]=(float*)malloc(3*sizeof(float));
        for(k=0;k<3;k++) rcm[i][k] = 0.0;
    }

    //get center of mass of triangle
    for(i=0;i<nMol;i++){
        float dx[3];
        for(k=0;k<3;k++) dx[k] = rru0[i*NN*3+k]-rr0[i*NN*3+k]; //use first bead as reference;
        float rTemp[3][3]; //consider only triangle, hence rTemp[3][3] and nor rTemp[NN][3];
        for(k=0;k<3;k++) rTemp[0][k] = rr0[i*NN*3+k];
        for(j=1;j<3;j++){
            for(k=0;k<3;k++) rTemp[j][k] = rru0[i*NN*3+3*j+k]-dx[k];
        }
        for(j=0;j<3;j++){
            for(k=0;k<3;k++) rcm[i][k]+= rTemp[j][k];
        }
        for(k=0;k<3;k++) rcm[i][k]/=3.0;
        
        //correct for periodicity
        for(k=0;k<3;k++){
            if(rcm[i][k]>pbcl[k]){
                rcm[i][k]-=2.0*pbcl[k];
            }
            else if (rcm[i][k] < -pbcl[k]){
                rcm[i][k]+=2.0*pbcl[k];
            }
        }
        
        //printf("%f %f %f\n",rcm[i][0],rcm[i][1],rcm[i][2]);
        
    }

    float factor = 4.0;
    int nBins = (int) factor * roundf(1.25*pbcl[0]);
    float dBins = 1.25*pbcl[0]/nBins/factor;
    float *rDist;
    rDist = (float*) malloc (nBins * sizeof(float));
    for(i=0;i<nBins;i++) rDist[i] = 0.0;
    int rTot=0;

    for(i=0;i<nMol;i++){
        for(j=i+1;j<nMol;j++){
            float dr2=0;
            for(k=0;k<3;k++){
                float dr;
                dr = fabsf(rcm[i][k] -rcm[j][k]);
                if (dr>pbcl[k]) dr = 2.0*pbcl[k] - dr;
                dr2+= dr*dr;
            }
            dr2 = sqrtf(dr2);
            //printf("%i %i %f\n",i,j,dr2);
            int ii = (int) floorf(dr2/dBins);
            if(ii<nBins) rDist[ii]++;
        }
    }
    
    FILE *fpw0;
    fpw0=fopen("gofr.txt","w");
    float rho_bulk = nMol/8.0/pbcl[0]/pbcl[1]/pbcl[2];
    for(i=0;i<nBins;i++){
        float dVol = 4.0*M_PI * (powf(dBins*(i+1.0),3.0)- powf(dBins*i,3.0))/3.0;
        fprintf(fpw0,"%f %f\n", dBins*(i+0.5),2.0*rDist[i]/dVol/rho_bulk/nMol);
    }
    fclose(fpw0);

}

