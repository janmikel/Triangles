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

#include "RgTensor.h"

int NN = 6; //number of beads in one molecule
float cluster_cutoff = 2.5; //based on g(r) where a peak of 1.47 is found
int min_cluster_size = 3; //based on inflection point of cluster size distribution

//Gives the cross product of two vectors
void CrossProduct(float *A, float *B, float *C){
//AXB=C

    C[0]= A[1]*B[2] - A[2]*B[1];
    C[1]= A[2]*B[0] - A[0]*B[2];
    C[2]= A[0]*B[1] - A[1]*B[0];
}

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

    int *cluster_index;
    cluster_index= (int*) malloc (nMol * sizeof(int));
    for(i=0;i<nMol;i++) cluster_index[i] = i;
    
    //loop through all contacts
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
            
            if(dr2<cluster_cutoff){
                if(cluster_index[i] < cluster_index[j]) cluster_index[j] = cluster_index[i]; //normal progresion
                else if ((cluster_index[i] > cluster_index[j])){
                    //revisit contacts and change indexes to new index
                    int old_index = cluster_index[i];
                    for(k=0;k<j;k++){
                        if(cluster_index[k] == old_index) cluster_index[k] = cluster_index[j];
                    }
                }
            }
        }
    }
    
    //get cluster size by counting cluster indexes
    int *cluster_size;
    cluster_size= (int*) malloc (nMol * sizeof(int));
    for(i=0;i<nMol;i++) cluster_size[i] = 0;
    
    for(i=0;i<nMol;i++){
        int index = cluster_index[i];
        cluster_size[index] ++;
    }
    
    FILE *fpw0;
    fpw0=fopen("orient.txt","w");
   //loop through clusters with size >= min_cluster_size 
   for(i=0;i<nMol;i++){ //loop through cluster indexes
        if(cluster_size[i]>= min_cluster_size){
            
            //initialize S tensor
            double **nSS,**SS; //uniaxial S tensor
            int ii,jj;
            nSS=matrix(1,3,1,3);
            SS=matrix(1,3,1,3);
            
            for(ii=0;ii<3;ii++){
                for(jj=0;jj<3;jj++){
                    nSS[ii+1][jj+1] = 0.0;
                    SS[ii+1][jj+1] = 0.0;
                }
            }
            
            //loop through all center of mass and determine molecules with the same cluster_index
            for(j=0;j<nMol;j++){
                if(i==cluster_index[j]){
                    //printf("%f %f %f %i\n",rcm[j][0],rcm[j][1],rcm[j][2],cluster_index[j]);
                    //get unit vector of rcm[j]
                    float A[3],B[3],C[3],mag;
                    for(k=0;k<3;k++){
                        //NN=1 and NN=0
                        A[k] = rru0[j*NN*3+3+k]-rru0[j*NN*3+k];
                        //NN=2 and NN=1
                        B[k] = rru0[j*NN*3+6+k]-rru0[j*NN*3+3+k];
                    }
                    CrossProduct(A,B,C);
                    
                    //get unit vector
                    mag=0.0;
                    for(k=0;k<3;k++){
                        mag+=C[k]*C[k];
                    }
                    mag=sqrtf(mag);
                    for(k=0;k<3;k++){
                        C[k]/=mag;
                    }
                    
                    //populate S tensor
                    for(ii=0;ii<3;ii++){
                        for(jj=0;jj<3;jj++){
                            double kdel = 0.0;
                            if(ii==jj) kdel = 1.0;
                            SS[ii+1][jj+1] += (C[ii]*C[jj] - kdel/3.0);
                            nSS[ii+1][jj+1]++;
                        }
                    }
                    
                }
            }
            
            //normalize S tensor
             for(ii=0;ii<3;ii++){
                for(jj=0;jj<3;jj++){
                    SS[ii+1][jj+1] /= nSS[ii+1][jj+1];
                }
            }
            
            
            //printf("Nss %i %i\n",(int) nSS[1][1], (int) nSS[2][3]);
            //get the Eigenvalues of the tensor
            double *d, **v;
            d = vectort(1,3);
            v= matrix(1,3,1,3);
            int nrot;
            jacobi(SS,3,d,v,nrot);
            eigsrt(d,v,3);
            
            fprintf(fpw0,"%i %f\n", (int) nSS[1][1],1.5*d[1]);
            
            free_vector(d,1,3);
            free_matrix(SS,1,3,1,3);
            free_matrix(nSS,1,3,1,3);
            free_matrix(v,1,3,1,3);
            
            
        }
    }
    fclose(fpw0);

}

