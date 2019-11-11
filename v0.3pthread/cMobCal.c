#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <pthread.h>

// Here I set a few constants
#ifndef  TEST
    #define  TEST   0
#endif

#ifndef  TEST_CTYPES
    #define  TEST_CTYPES  0
#endif

#define  N_CORES    2

#define  CANG   (180./M_PI)
#define  XE     1.60217733e-19
#define  XK     1.380658e-23
#define  XN     6.0221367e23
#define  XEO    8.854187817e-12
#define  XMV    0.02241410
#define  EO     (1.34e-03*(XE))
#define  RO     3.043e-10
#define  RO2    ((RO)*(RO))

#if TEST
    // INP should be at least 5
    #define  IPR    2
    #define  INP    5   
    #define  ITN    3 
    #define  IMP    4
    #define  INUM   30
    #define  INOR   20
#else
    #define  IPR    1000
    #define  INP    40   
    #define  ITN    10 
    #define  IMP    25
    #define  INUM   250000
    #define  INOR   30

#endif

#define  CMIN   0.0005
#define  SW1    0.00005
#define  SW2    0.005 
#define  DTSF1  0.5
#define  DTSF2  0.1
#define  INWR   1
#define  IFAIL  100
#define  NTHETA 33.3
#define  NPHI   64.250
#define  NGAMMA       134.30
#define  TRAJ_LOST    30000        
#define  TEMPERATURE  298
#define  IBST_MAX     500
#define  MAX_FILENAME 30
#define  MAX_LINE     100
#define  XL_INIT      1e+6

// Here I set some constant for printing verbosity.
// DIFF: Most of the following switches have been implemented as macros, to improve perfomance
// NOTE: setting PRINT_IT has actually *no* effect in the original mobcal (broken?). In this new version, it works as expected
//       When comparing to original mobcal, keep PRINT_IT=0
#define PRINT_IT  0
#define PRINT_IP  0
#define PRINT_IGS 0
#define PRINT_IU1 0
#define PRINT_IU2 0
#define PRINT_IU3 0


// DIFF: The following PRINT switches have been removed from this list as the latest
// version of mobcal set them in the main function as needed. Thus the value set here
// has little to none effect
//#define PRINT_IM2 0
//#define PRINT_IM4 0

struct Molecule{
    /* 
    * Structure containing information about the system under investigation
    * 
    * Attributes
    * ----------
    * TotalMass : total mass of the molecule (excluding the gas particle)
    * dipole : dipole of the molecule (excluding the gas particle)
    * totalCharge : total charge of the molecule (excluding the gas particle)
    * totalAbsoluteCharge : sum of the absolute value of the charges (excluding the gas particle)
    * mGas : mass of the gas particle
    * mu : reduced mass of the system (molecule + gas particle)
    * invMu: 1/mu. divisions is an expensive operation, so we compute it once here and reuse
    * mobility : mobility coefficient of the molecule (excluding the gas particle)
    * temperature : temperature of the system in Kelvin
    * 
    */    

    double TotalMass;
    double dipole;
    double totalCharge;
    double totalAbsoluteCharge;
    double mGas;
    double mu;     
    double invMu;  
    double mobility;
    double temperature;  
};

struct Molecule InitMolecule(void){
    /*
    * This function allocs and initialize the Molecule structure that contains informations
    * about the system investigated
    * 
    * Returns
    * -------
    * struct Molecules: newly create struct with data. 
    * See definition of struct for the meaning of each attribute.
    */
    
    struct Molecule mol;

    mol.TotalMass=0;
    mol.dipole=0;
    mol.totalCharge=0;
    mol.totalAbsoluteCharge=0;
    mol.mGas=0;
    mol.mu=0;     
    mol.invMu=0;  
    mol.mobility=0;
    mol.temperature=TEMPERATURE;  

    return mol;
}


struct LennardJones{
    /* 
    * Structure containing information for calculating Lennard Jones interactions
    * The code defines the potential as U = 4 * epsilon *[ (sigma/r)^12 - (sigma/r)^6 ]
    * 
    * Attributes
    * ----------
    * eolj : list of epsilon (one for each atom)
    * rolj : list of sigmas (one for each atom)
    * n : number of atoms
    * romax : largest value in rolj
    * 
    * The following attributes are derived from the previous.
    * As they are used over and over, we calculate them once and reuse.
    * 
    * eox4 : 4 * eolj (one for each atom)
    * ro2lj : square of rolj (one for each atom)
    * ro6lj : rolj^6 (one for each atom)
    * ro12lj : rolj^12 (one for each atom)
    * dro6 : 6 * ro6lj (one for each atom)
    * dro12 : 12 * ro12lj (one for each atom)    
    */    
    double * eolj;
    double * rolj;
    double *eox4;
    double *ro2lj;
    double *ro6lj;
    double *ro12lj;
    double *dro6;
    double *dro12;
    double romax; 
    size_t n;
};

double* AllocDoubleVec(size_t n){
    /*
    * This function allocs the memory for a vector
    * 
    * Arguments
    * ---------
    * n: the length of the array
    * 
    * Returns
    * -------
    * double *: a pointer to a vector
    */

    double *v=NULL;

    v=(double*) malloc(n*sizeof(double));
    if (v==NULL){
        puts("I could not create the vector");
        exit(EXIT_FAILURE);
    }
    return v;
}

unsigned int* AllocIntVec(size_t n){
    /*
    * This function allocs the memory for a vector
    * 
    * Arguments
    * ---------
    * n: the length of the array
    * 
    * Returns
    * -------
    * double *: a pointer to a vector
    */

    unsigned int *v=NULL;

    v=(unsigned int*) malloc(n*sizeof(int));
    if (v==NULL){
        puts("I could not create the vector");
        exit(EXIT_FAILURE);
    }
    return v;
}

double** AllocDoubleMat(size_t rows,size_t columns){
    /*
    * This function allocs the memory for a matrix
    * 
    * Arguments
    * ---------
    * rows: the number of rows
    * columns: the number of columns
    * 
    * Returns
    * -------
    * double **: a pointer to a matrix 
    */

    double **v=NULL;
    double *array=NULL;
    unsigned int i;

    v=(double**) malloc(rows*sizeof(double*));
    array=(double*) malloc(rows*columns*sizeof(double));
    if (v==NULL || array==NULL){
        puts("I could not create the matrix");
        exit(EXIT_FAILURE);
    }

    for(i=0;i<rows;i++) v[i]=&array[i*columns];
        
    return v;
}

void FreeDoubleMat(double **v){
    /*
    * This function frees the memory allocated for a matrix
    * 
    * Arguments
    * ---------
    * v : pointer to the matrix
    * 
    */
    
    free(v[0]);
    free(v);

}


struct HardSphere {
    /* 
    * Structure containing information for calculating hard sphere scattering
    * 
    * Attributes
    * ----------
    * rhs : list of radii  (one for each atom)
    * rhs2 : square of rhs. It is used over and over, so we compute and reuse to speed up calculations
    * 
    */    

    double *rhs;
    double *rhs2;

};

void FreeHS(struct HardSphere *hs){
    /*
    * Function that free the memory allocated for this structure
    * 
    * Arguments
    * ---------
    * hs : structure to be freed
    *
    */

    free(hs->rhs);
    free(hs->rhs2);
}

struct LennardJones SetLJ(double* eolj,double* rolj, size_t n){
    /*
    * This function allocs the memory for storing Lennard-Jones parameters
    * The code defines the potential as U = 4 * epsilon *[ (sigma/r)^12 - (sigma/r)^6 ]
    * 
    * Arguments
    * ---------
    * eolj: vector with the list of epsilon (one per atom)
    * rolj: vector with the list of sigma (one per atom)
    * n: number of atoms
    * NOTE: n *must* be the length of both eolj and rolj. Thus, eolj and rolj have same length.
    * function does *not* check for that.
    * 
    * Returns
    * -------
    * struct LennardJones: newly create struct with data. 
    * See definition of struct for the meaning of each attribute.
    */
    
    struct LennardJones lj;
    unsigned int i;

    lj.n=n;
    lj.eolj=eolj;
    lj.rolj=rolj;
    
    lj.eox4=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.eox4[i]=4.*lj.eolj[i];
    }

    lj.ro2lj=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.ro2lj[i]=lj.rolj[i]*lj.rolj[i];
    }

    lj.ro6lj=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.ro6lj[i]=lj.ro2lj[i]*lj.ro2lj[i]*lj.ro2lj[i];
    }

    lj.ro12lj=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.ro12lj[i]=lj.ro6lj[i]*lj.ro6lj[i];
    }

    lj.dro6=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.dro6[i]=6.*lj.ro6lj[i];
    }

    lj.dro12=AllocDoubleVec(n);
    for(i=0;i<n;i++){
        lj.dro12[i]=12.*lj.ro12lj[i];
    }

    lj.romax=lj.rolj[0];
    for(i=1;i<n;i++){
        if (lj.rolj[i]>lj.romax) lj.romax=lj.rolj[i];
    }

    return lj;
}

void FreeLj(struct LennardJones* lj){
    /*
    * Function that free the memory allocated for this structure
    * 
    * Arguments
    * ---------
    * lj: structure to be freed
    *
    */
    free(lj->eolj);
    free(lj->rolj);
    free(lj->eox4);
    free(lj->ro2lj);
    free(lj->ro6lj);
    free(lj->ro12lj);
    free(lj->dro6);
    free(lj->dro12);
}

struct atom{
    /* 
    * Structure containing information about a specific atomtype
    * 
    * Attributes
    * ----------
    * name : name of the atom
    * mass : atomic mass
    * epsilonLj,sigmaLJ : Lennard-Jones parameters
    * rhs : parameters for hard sphere scattering
    * id : atomic number
    * 
    */    

    char name[25];
    double mass;
    double epsilonLJ;
    double sigmaLJ;
    double rhs;
    unsigned int id;
};

struct arg_mc{
    double **coord;
    struct HardSphere* hs;
    #if TEST_CTYPES
    double (*RanNum)(void);
    #else
    double (*RanNum)(const gsl_rng *r);
    #endif
    const gsl_rng *r;
    double *crb;
    unsigned int * imm;
    double *crof;
    size_t n_atoms;    
};

#if TEST_CTYPES
void SetArgmc(struct arg_mc *args,double **coord,struct HardSphere* hs,double (*RanNum)(void),double *crb,unsigned int * imm,double *crof,size_t n_atoms){
args->r=NULL;
#else
void SetArgmc(struct arg_mc *args,double **coord,struct HardSphere* hs,double (*RanNum)(const gsl_rng *r),const gsl_rng *r,double *crb,unsigned int * imm,double *crof,size_t n_atoms){
args->r=r;
#endif
args->coord=coord;
args->hs=hs;
args->RanNum=RanNum;
args->crb=crb;
args->imm=imm;
args->crof=crof;
args->n_atoms=n_atoms;
}

struct arg_MD{
    double **coord;
    const struct LennardJones *lj;
    #if TEST_CTYPES
    double (*RanNum)(void);
    #else
    double (*RanNum)(const gsl_rng *r);
    #endif
    const gsl_rng *r;
    const double *charge;
    const struct Molecule *mol;
    size_t n_atoms;    
    unsigned int t_id;
    double *b2max;
    double *vList;
    int  *ifailc;
    double tst;
    double *wgst;
    double *pgst;
    double *om11st;
    double *om12st;
    double *om13st;
    double *om22st;
    double *q1st;
    double *q2st;
};

#if TEST_CTYPES
void SetArgMD(struct arg_MD *args,unsigned int t_id,double **coord,const struct LennardJones* const lj,const double* const charge,double (*RanNum)(void),size_t n_atoms,double *b2max, double *vList,int *ifailc,double tst,double *om11st,double *om12st,double *om13st,double *om22st,double *q1st,double *q2st){
args->r=NULL;
#else
void SetArgMD(struct arg_MD *args,unsigned int t_id,double **coord,const struct LennardJones* const lj,const double* const charge,const struct Molecule *const mol,double (*RanNum)(const gsl_rng *r),const gsl_rng *r,size_t n_atoms,double *b2max, double *vList,double *wgst,double *pgst,int *ifailc,double tst,double *om11st,double *om12st,double *om13st,double *om22st,double *q1st,double *q2st){
args->r=r;
#endif

args->coord=coord;
args->lj=lj;
args->RanNum=RanNum;
args->n_atoms=n_atoms;
args->charge=charge;
args->t_id=t_id;
args->b2max=b2max;
args->ifailc=ifailc;
args->vList=vList;
args->tst=tst;
args->wgst=wgst;
args->pgst=pgst;
args->om11st=om11st;
args->om12st=om12st;
args->om13st=om13st;
args->om22st=om22st;
args->q1st=q1st;
args->q2st=q2st;
args->mol=mol;
}

const struct atom AtomList[]={
    /*
    * Here is a list with the atom recognized by MobCal.
    * For new atoms, add a line here.
    * The software defaults to carbon (and issues a warning) if the
    * atomtype specified is not found here.
    */
    {"Hydrogen",1.008,0.65e-03*XE,2.38e-10,2.2e-10,1},
    {"Carbon",12.01,1.34e-3*XE,3.043e-10,2.7e-10,12},
    {"Nitrogen",14.01,1.34e-3*XE,3.043e-10,2.7e-10,14},
    {"Oxygen",16.00,1.34e-3*XE,3.043e-10,2.7e-10,16},
    {"Sodium",22.99,0.0278e-3*XE,(3.97/1.12246)*1e-10,2.853e-10,23},
    {"Silicon",28.09,1.35e-3*XE,3.5e-10,2.95e-10,28},
    {"Sulfur",32.06,1.35e-3*XE,3.5e-10,3.5e-10,32},
    {"Iron",55.85,1.35e-3*XE,3.5e-10,3.5e-10,56}
};

#define ATOM_LIST_LENGTH (sizeof(AtomList)/sizeof(struct atom))

void AssignAtom(unsigned int id, double *mass,double *epsilonLJ,double *sigmaLJ, double *rhs){
    /*
    * This function assign atom properties based on atom atomic number
    * This function comapares the mass found in teh file with the database and set all
    * the internal variable for that atom.
    * If the atom is NOT in the database it is defaulted to a carbon atom and a warning is issued.
    * Whether using a carbon atom as default is a good idea or not, I will let you decide
    * 
    * Arguments
    * ---------
    * id : atomic number
    * 
    * Returns
    * -------
    * mass : atom mass
    * epsilonLJ, sigmaLJ : Lennard-Jones parameters
    * rhs : hard shpere radius 
    */    

    unsigned int i;

    for(i=0;i<ATOM_LIST_LENGTH;i++){
        if(AtomList[i].id==id){
            *mass=AtomList[i].mass;
            *epsilonLJ=AtomList[i].epsilonLJ;
            *sigmaLJ=AtomList[i].sigmaLJ;
            *rhs=AtomList[i].rhs;
            return;
        }
    }
    // if the atomtype does not match
    puts("****** WARNING   WARNING   WARNING ******");
    puts("One of the atoms listed is unknown, assuming carbon");
    puts("****** WARNING   WARNING   WARNING ******");

    for(i=0;i<ATOM_LIST_LENGTH;i++){
        if(AtomList[i].id==12){
            *mass=AtomList[i].mass;
            *epsilonLJ=AtomList[i].epsilonLJ;
            *sigmaLJ=AtomList[i].sigmaLJ;
            *rhs=AtomList[i].rhs;
            return;
        }
    }

}

struct inputFileData{
    /*
    * This structure contains some general information obtained from the
    * input file. Ideally, no matter what the actually file format looks like,
    * the function reading it should populate thise structure.
    * 
    * Attributes
    * ----------
    * label : a user-provided label for the simulation. Limit to MAX_LINE
    * n_atoms : number of atoms per molecules
    * n_coords : Sets of coordinates. Ideally, different structures of the same molecules. 
    *            *not* different orientation as the program takes care of that.
    * units : the units used for the coordinates of the atoms. AU=atomic units, ANG=angstroms
    * chargeUnits : How the charges are generated. NONE = no charges, EQUAL = all charges equal 1/n_Atoms, 
    *               CALC = read from file (an extra column must be present in the file)
    */

    char label[MAX_LINE];
    unsigned int n_atoms;
    unsigned int n_coord;
    enum units{AU,ANG} unit;
    enum chargeUnits{NONE,EQUAL,CALC} cUnit;
    double cFact;
};

struct AtomData{
    /*
    * This structure contains about each atom in the molecules
    * 
    * Attributes
    * ----------
    * coord : coordinates
    * charge : charges
    * mass : masses
    */

    double **coord;
    double *charge;
    double *mass;
};

void FreeAtomData(struct AtomData* atoms){
    /*
    * Function that free the memory allocated for this structure
    * 
    * Arguments
    * ---------
    * atoms : structure to be freed
    *
    */
   
    FreeDoubleMat(atoms->coord);
    free(atoms->charge);
    free(atoms->mass);
}


void RotZ(const double *const in,double *out,double angle){
    double Rz[3][3];
    /*
    * This function rotates a vector around the Z axis.
    * 
    * Arguments
    * ---------
    * in : vector to be rotated
    * angle : angle of the rotation
    * 
    * Returns
    * -------
    * out : rotated vector
    */    

    // Define the rotation matrix
    Rz[0][0]=cos(angle); Rz[0][1]=-sin(angle); Rz[0][2]=0;
    Rz[1][0]=sin(angle); Rz[1][1]=cos(angle);  Rz[1][2]=0;
    Rz[2][0]=0;          Rz[2][1]=0;           Rz[2][2]=1;

    // Apply the rotation
    out[0]=in[0]*Rz[0][0]+in[1]*Rz[0][1]+in[2]*Rz[0][2];
    out[1]=in[0]*Rz[1][0]+in[1]*Rz[1][1]+in[2]*Rz[1][2];
    out[2]=in[0]*Rz[2][0]+in[1]*Rz[2][1]+in[2]*Rz[2][2];

}

void RotX(const double *const in,double *out,double angle){
    double Rx[3][3];
    /*
    * This function rotates a vector around the X axis.
    * 
    * Arguments
    * ---------
    * in : vector to be rotated
    * angle : angle of the rotation
    * 
    * Returns
    * -------
    * out : rotated vector
    */    

    // Define the rotation matrix
    Rx[0][0]=1;          Rx[0][1]=0;           Rx[0][2]=0;
    Rx[1][0]=0;          Rx[1][1]=cos(angle);  Rx[1][2]=-sin(angle);
    Rx[2][0]=0;          Rx[2][1]=sin(angle);  Rx[2][2]=cos(angle);

    // Apply the rotation
    out[0]=in[0]*Rx[0][0]+in[1]*Rx[0][1]+in[2]*Rx[0][2];
    out[1]=in[0]*Rx[1][0]+in[1]*Rx[1][1]+in[2]*Rx[1][2];
    out[2]=in[0]*Rx[2][0]+in[1]*Rx[2][1]+in[2]*Rx[2][2];

}


void Rotate(double **InVector,double **OutVector,double theta,double phi,double gamma, size_t n_atoms){
    /*
    * This function rotates the molecule under investigation. 
    * The molecules is rotated around Z for theta radians,
    * then -phi radian around X and finally  gamma radians around Z
    * 
    * Arguments
    * ---------
    * InVector : coordinates of the molecules to be rotated
    * theta,phi,gamma : angles of the rotation
    * n_atoms : number of atoms within the molecule
    * fp : file for logging data
    * 
    * Returns
    * -------
    * OutVector : rotated molecule
    */    

    unsigned int i;
    double tmpVector[3];

    // Rotate vector
    for(i=0;i<n_atoms;i++){
        RotZ(InVector[i],OutVector[i],theta);
        RotX(OutVector[i],tmpVector,-phi);
        RotZ(tmpVector,OutVector[i],gamma);
    }
}

// DIFF : The original code allows new masses, but retains the same lj and rhs parameters. I read only positions and charges
void ncoord(struct inputFileData *inData,struct Molecule *const mol, struct AtomData *atoms, double *asymp, FILE *fpIn, FILE *fpOut){

    /*
    * Functions for updating the configuration analyzed from file.
    * 
    * Arguments
    * ---------
    * inData : input parameters read by fcoord()
    * mol : molecular properties
    * fpIn : file with the data to be read
    * fpOut: file for logging
    * 
    * Returns
    * -------
    * atoms : contains updated information about each atom in the molecule
    * 
    */    
    
    char buffer[MAX_LINE];
    char tmp[MAX_LINE];
    unsigned int i,j,k;

    #if (PRINT_IU2)
    unsigned int n;
    #endif

    double **rotCoords, cMass[3]={0.,0.,0.},theta,phi,gamma,xyzsum,yzsum,hold;

    //There is an unused line separating data sets. We read and trash
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: Missing new dataset ");
        exit(EXIT_FAILURE);
    }else if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)

    // Now let's read the data
    // Follows the information about atoms : x, y, z, mass, charge(optional)

    for(i=0;i<inData->n_atoms;i++){
        if (fgets(buffer,MAX_LINE,fpIn)==NULL){
            fclose(fpIn);
            fclose(fpOut);
            perror("Error: Missing line with coordinates ");
            exit(EXIT_FAILURE);
        }else{
            strcpy(tmp,strtok(buffer," \n"));
            atoms->coord[i][0]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            atoms->coord[i][1]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            atoms->coord[i][2]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            //Data Not Used
            if(inData->cUnit==CALC) {
                strcpy(tmp,strtok(NULL," \n"));
                atoms->charge[i]=strtod(tmp,NULL);
            //printf("%d %lf %lf %lf\n",i,coords[i][0],coords[i][1],coords[i][2]);    
            }
        }
    }
    if(inData->unit==AU) for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) atoms->coord[i][0]*=0.52917706;

    // The origin of the reference frame is the center of the molecules. We calculate it and shift the coordinates
    for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) cMass[j]+=atoms->mass[i]*atoms->coord[i][j];
    for (j=0;j<3;j++) cMass[j]/=mol->TotalMass;
    for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) atoms->coord[i][j]=(atoms->coord[i][j]-cMass[j])*1e-10*inData->cFact;

    // compute the structural asymmetry parameter
    *asymp=0.;
    theta=0;

    rotCoords=AllocDoubleMat(inData->n_atoms,3);
    for(i=0;i<=360;i+=2){
        for(j=0;j<=180;j+=2){
            xyzsum=0.;
            yzsum=0.;
            gamma=(double)i/CANG;
            phi=(double)j/CANG;
            Rotate(atoms->coord,rotCoords,theta,phi,gamma,inData->n_atoms);
            
            #if (PRINT_IU2 || PRINT_IU3)
            fprintf(fpOut,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
            #endif
            
            #if (PRINT_IU2)
            fprintf(fpOut,"         initial coordinates                        new coordinates\n\n");
            for(n=0;n<inData->n_atoms;n++) fprintf(fpOut," % -.4E % -.4E % -.4E      % -.4E % -.4E % -.4E\n",atoms->coord[n][0],atoms->coord[n][1],atoms->coord[n][2],rotCoords[n][0],rotCoords[n][1],rotCoords[n][2]);
            #endif

            for (k=0;k<inData->n_atoms;k++){
                xyzsum+=sqrt(rotCoords[k][0]*rotCoords[k][0]+rotCoords[k][1]*rotCoords[k][1]+rotCoords[k][2]*rotCoords[k][2]);
                yzsum+=sqrt(rotCoords[k][1]*rotCoords[k][1]+rotCoords[k][2]*rotCoords[k][2]);
            }
            hold=((M_PI_4)*xyzsum)/yzsum;
            if (hold>(*asymp)) *asymp=hold;
        }
    }

    FreeDoubleMat(rotCoords);
}

//return atoms,ion,InputData,asymp,lj,hs,fin,fout
void fcoord(const char *inFile,struct inputFileData *inData,struct LennardJones* lj,struct HardSphere* hs,struct Molecule *const mol, struct AtomData * const atoms, double *asymp, FILE *fpIn, FILE *fpOut){
    /*
    * Functions for updating the configuration analyzed from file.
    * 
    * Arguments
    * ---------
    * inFile : name of the input file
    * fpIn : file with the data to be read
    * fpOut : file for logging
    * 
    * Returns
    * -------
    * inData : several parameters needed for running the simulation
    * lj : lennard-jones parameters
    * hs : hards sphere parameters
    * mol : molecular properties
    * atoms : contains updated information about each atom in the molecule
    * asymp : structural asymmetry parameter
    * 
    */    

    char buffer[MAX_LINE];
    char tmp[MAX_LINE];
    unsigned int i,j,*massID,k;
    double **rotCoords,*charges,*masses,invCharge,*eolj,*rolj,*rhs,*rhs2,cMass[3]={0.,0.,0.},theta,phi,gamma,xyzsum,yzsum,hold;

    fprintf(fpOut," input file name = %s\n",inFile);
    
    // The standard MobCal files are position based. The next is the label used for the system
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: The label was not found ");
        exit(EXIT_FAILURE);
    }else {
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
        strcpy(inData->label,strtok(buffer,"\n"));  // remove trailing newlines
    }
    fprintf(fpOut," input file label = %s\n",inData->label);    

    //number of sets of coordinates
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: The number of coordinates was not found ");
        exit(EXIT_FAILURE);
    }else{        
        inData->n_coord=strtoul(buffer,NULL,10);
        //printf("%s => %d",buffer,inData->n_coord);
        if(!inData->n_coord){
        fclose(fpIn);
        fclose(fpOut);
        fprintf(stderr,"Error: The number of coordinates must be > 0\n");
        exit(EXIT_FAILURE);            
        }else fprintf(fpOut," number of coordinate sets =%5d\n",inData->n_coord);       
    }
    if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
    
    // number of atoms in the molecule
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: The number of atoms was not found ");
        exit(EXIT_FAILURE);
    }else{
        inData->n_atoms=strtoul(buffer,NULL,10);
        if(!inData->n_atoms){
        fclose(fpIn);
        fclose(fpOut);
        fprintf(stderr,"Error: The number of atoms must be > 0\n");
        exit(EXIT_FAILURE);            
        }else fprintf(fpOut," number of atoms =%4d\n",inData->n_atoms);       
    }
    if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)

    // What units are used?
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: the units was not found ");
        exit(EXIT_FAILURE);
    }else{
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)        
        strcpy(tmp,strtok(buffer," \n")); //remove trailing new lines 
        if(!strcmp(tmp,"au")) inData->unit=AU;
        else if (!strcmp(tmp,"ang")) inData->unit=ANG;
        else {
        fclose(fpIn);
        fclose(fpOut);
        fprintf(stderr,"Error: Type of coordinates not supported\n");
        exit(EXIT_FAILURE);            
        }        
    }

    switch(inData->unit){
        case 0: fprintf(fpOut," coordinates in atomic units\n");
            break;       
        case 1: fprintf(fpOut," coordinates in angstroms\n");
            break;       
    }
    
    // Here we write the correction factor *before* charges.
    // Not sure why, but this is done in the original code
    // and I am retaining compatibility

    // Lets read charge units
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: the units for the charges was not found ");
        exit(EXIT_FAILURE);
    }else{
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
        strcpy(tmp,strtok(buffer," \n")); //remove trailing new lines 
        if(!strcmp(tmp,"equal")) inData->cUnit=EQUAL;
        else if (!strcmp(tmp,"calc")) inData->cUnit=CALC;
        else if (!strcmp(tmp,"none")) inData->cUnit=NONE;
        else {
        fclose(fpIn);
        fclose(fpOut);
        fprintf(stderr,"Error: Type of charges not supported\n");
        exit(EXIT_FAILURE);            
        }        
    }

    //Let's read the correction factor;
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        fclose(fpOut);
        perror("Error: The correction factor was not found ");
        exit(EXIT_FAILURE);
    }else{
        inData->cFact=strtod(buffer,NULL);
        if(!inData->cFact){
        fclose(fpIn);
        fclose(fpOut);
        fprintf(stderr,"Error: The correction factor must be > 0\n");
        exit(EXIT_FAILURE);            
        }else fprintf(fpOut," correction factor for coordinates =% -.4E\n",inData->cFact);
    }
    if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)

    // Now we can print information about the charges
    switch(inData->cUnit){
        case 0: fprintf(fpOut," using no charge - only LJ interactions\n");
            break;       
        case 1: fprintf(fpOut," using a uniform charge distribution\n");
            break;       
        case 2: fprintf(fpOut," using a calculated (non-uniform) charge distribution\n");
            break;       
    }

    // Follow the information about atoms : x, y, z, mass, charge(optional)
    atoms->coord=AllocDoubleMat(inData->n_atoms,3);
    massID=AllocIntVec(inData->n_atoms);
    masses=AllocDoubleVec(inData->n_atoms);
    charges=AllocDoubleVec(inData->n_atoms);
    eolj=AllocDoubleVec(inData->n_atoms);
    rolj=AllocDoubleVec(inData->n_atoms);
    rhs=AllocDoubleVec(inData->n_atoms);
    rhs2=AllocDoubleVec(inData->n_atoms);

    for(i=0;i<inData->n_atoms;i++){
        if (fgets(buffer,MAX_LINE,fpIn)==NULL){
            fclose(fpIn);
            fclose(fpOut);
            perror("Error: Missing line with coordinates ");
            exit(EXIT_FAILURE);
        }else{
            //printf("%d\n",i);
            strcpy(tmp,strtok(buffer," \n"));
            atoms->coord[i][0]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            atoms->coord[i][1]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            atoms->coord[i][2]=strtod(tmp,NULL);
            strcpy(tmp,strtok(NULL," \n"));
            massID[i]=(int)round(strtod(tmp,NULL));
            AssignAtom(massID[i], &masses[i], &eolj[i], &rolj[i], &rhs[i]);
            rhs2[i]=rhs[i]*rhs[i];            
            if(inData->cUnit==CALC) {
                strcpy(tmp,strtok(NULL," \n"));
                charges[i]=strtod(tmp,NULL);
            }        
        }
    }

    if(inData->cUnit==EQUAL){
        invCharge=1./(double) inData->n_atoms;
        for(i=0;i<inData->n_atoms;i++) charges[i]=invCharge;
    }
    
    if(inData->cUnit==NONE) for(i=0;i<inData->n_atoms;i++) charges[i]=0.;

    if(inData->unit==AU) for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) atoms->coord[i][0]*=0.52917706;

    //And we finally link everything to its structure 
    //atoms->coord=coords;
    atoms->mass=masses;
    atoms->charge=charges;
    *lj=SetLJ(eolj,rolj,inData->n_atoms);
    hs->rhs=rhs;
    hs->rhs2=rhs2;

    // Now we compute some global properties of the entire molecules in addition to atom by atom properties just read

    //total charge:
    mol->totalCharge=0;
    for (i=0;i<inData->n_atoms;i++) mol->totalCharge+=atoms->charge[i];

    //total absolute charge:
    mol->totalAbsoluteCharge=0.;
    for (i=0;i<inData->n_atoms;i++)  mol->totalAbsoluteCharge+=fabs(atoms->charge[i]);

    fprintf(fpOut," total charge =% -.4E\n total absolute charge =% -.4E\n",mol->totalCharge,mol->totalAbsoluteCharge);

    //total mass:
    mol->TotalMass=0.;
    for (i=0;i<inData->n_atoms;i++)  mol->TotalMass+=atoms->mass[i];

    // DIFF : Original MobCal only prints the charges is model==calc (see next 2 lines). I prefer always: helps catching bugs
    // DIFF : fortran uses D instead of E to point out that you are usingthe scientific notation to represent an integer number and not a float
    fprintf(fpOut," mass of ion =% -.4E\n",mol->TotalMass);

    #if PRINT_IU1 
    fprintf(fpOut,"\n         initial coordinates         mass   charge         LJ parameters\n\n");
    #endif

    // The origin of the reference frame is the center of the molecules. We calculate it .... 
    for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) cMass[j]+=atoms->mass[i]*atoms->coord[i][j];
    for(j=0;j<3;j++) cMass[j]/=mol->TotalMass;
    fprintf(fpOut," center of mass coordinates = % -.4E,% -.4E,% -.4E\n",cMass[0],cMass[1],cMass[2]);
    // ... and shift the coordinates directly into the final array
    for (i=0;i<inData->n_atoms;i++) for(j=0;j<3;j++) atoms->coord[i][j]=(atoms->coord[i][j]-cMass[j])*1e-10*inData->cFact;

    #if PRINT_IU1 
    for (i=0;i<inData->n_atoms;i++) fprintf(fpOut," % -.4E % -.4E % -.4E %3d % -.4E % -.4E % -.4E\n",atoms->coord[i][0],atoms->coord[i][1],atoms->coord[i][2],massID[i],atoms->charge[i],lj->eolj[i]/XE,lj->rolj[i]*1e10);
    fprintf(fpOut,"\n\n");
    #endif

    free(massID);
    
    if(inData->n_coord==1) fclose(fpIn);

    // compute the structural asymmetry parameter
    *asymp=0.;
    theta=0;

    rotCoords=AllocDoubleMat(inData->n_atoms,3);
    for(i=0;i<=360;i+=2){
        for(j=0;j<=180;j+=2){
            xyzsum=0.;
            yzsum=0.;
            gamma=(double)i/CANG;
            phi=(double)j/CANG;
            Rotate(atoms->coord,rotCoords,theta,phi,gamma,inData->n_atoms);

            #if (PRINT_IU2 || PRINT_IU3)
            fprintf(fpOut,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
            #endif
            
            #if (PRINT_IU2)
            fprintf(fpOut,"         initial coordinates                        new coordinates\n\n");
            for(n=0;n<inData->n_atoms;n++) fprintf(fpOut," % -.4E % -.4E % -.4E      % -.4E % -.4E % -.4E\n",atoms->coord[n][0],atoms->coord[n][1],atoms->coord[n][2],rotCoords[n][0],rotCoords[n][1],rotCoords[n][2]);
            #endif

            for (k=0;k<inData->n_atoms;k++){
                xyzsum+=sqrt(rotCoords[k][0]*rotCoords[k][0]+rotCoords[k][1]*rotCoords[k][1]+rotCoords[k][2]*rotCoords[k][2]);
                yzsum+=sqrt(rotCoords[k][1]*rotCoords[k][1]+rotCoords[k][2]*rotCoords[k][2]);
            }
            hold=((M_PI_4)*xyzsum)/yzsum;
            if (hold>(*asymp)) *asymp=hold;
        }
    }

    FreeDoubleMat(rotCoords);
}

void dljpot(const double* const position, double ** coord, const struct LennardJones* const lj,const double* const charge, double dipol, size_t n_atoms,double totCharge,double *pot, double dpot[3],double *dmax){
    /*
    * Function that computes the Lennard-Jones potential + dipolar interactions
    * The code defines the potential as U = 4 * epsilon *[ (sigma/r)^12 - (sigma/r)^6 ]
    * 
    * Arguments
    * ---------
    * position : vector with the position of the gas particle. *must* be of dimesion = 3
    * coord : matrix with the position of each atom. Each row is a different atom an contains the coordinates x,y,z
    * lj : structure with info about the Lennard-Jones parameters
    * charge : vector with the charges of each atom in the molecules
    * dipol : total dipole of the molecule
    * n_atoms : number of atoms within the molecule
    * totCharge : total charge of the molecule
    * 
    * Returns
    * -------
    * pot : total potential energy (LJ + dipole)
    * dpot : vector of dimension 3. derivative of the potential along each direction x,y,z. NOTE: the force is -dpot (*minus* dpot)
    * dmax : minimum distance between the particle of the gas and the molecule (namely, closest atom of the molecule)
    *    
    */

    double e00=0., de00x=0., de00y=0., de00z=0.,rx=0.,ry=0.,rz=0.,sum1=0.,sum2=0.,sum3=0.,sum4=0.,sum5=0.,sum6=0.;
    double data[3],data2[3],r,inv_r,inv_r2,r3,r5,r6,r8,r12,r14,qr3,qr5,de00;

    unsigned int i,j;

    *pot=0.;
    *dmax=lj->romax;
    for(i=0;i<n_atoms;i++){
        for(j=0;j<3;j++){
            data[j]=position[j]-coord[i][j];
            data2[j]=data[j]*data[j];    
        }
        r=sqrt(data2[0]+data2[1]+data2[2]);
        if (r<*dmax) *dmax=r;
        inv_r=1./r;
        inv_r2=inv_r*inv_r;
        r3=inv_r2*inv_r;
        r5=inv_r2*r3;
        r6=r5*inv_r;
        r8=r5*r3;
        r12=r6*r6;
        r14=r12*inv_r2;
        
        e00=e00+lj->eox4[i]*((lj->ro12lj[i]*r12)-(lj->ro6lj[i]*r6));
        de00=lj->eox4[i]*((lj->dro6[i]*r8)-(lj->dro12[i]*r14)); //no need to store de00? use and forget?
        de00x=de00x+de00*data[0];
        de00y=de00y+de00*data[1];
        de00z=de00z+de00*data[2];
        if (totCharge){
            qr3=charge[i]*r3;
            qr5=-3.*charge[i]*r5;
            rx=rx+qr3*data[0];
            ry=ry+qr3*data[1];
            rz=rz+qr3*data[2];

            sum1=sum1+qr3+data2[0]*qr5;
            sum2=sum2+data[0]*data[1]*qr5;
            sum3=sum3+data[0]*data[2]*qr5;
            sum4=sum4+qr3+data2[1]*qr5;
            sum5=sum5+data[1]*data[2]*qr5;
            sum6=sum6+qr3+data2[2]*qr5;
        }

    }
    // Total potential energy
    *pot=e00-(dipol*((rx*rx)+(ry*ry)+(rz*rz)));
    //  Derivatives of the potential with respect coords = derivatives of the  momenta  with respect time
    dpot[0]=de00x-(dipol*((2.0*rx*sum1)+(2.*ry*sum2)+(2.*rz*sum3)));
    dpot[1]=de00y-(dipol*((2.0*rx*sum2)+(2.*ry*sum4)+(2.*rz*sum5)));
    dpot[2]=de00z-(dipol*((2.0*rx*sum3)+(2.*ry*sum5)+(2.*rz*sum6)));

}

void deriv(const double* const position, const double* const momentum, double ** coord,const double invMass,const struct LennardJones* const lj,const double* const charge, double dipol, size_t n_atoms,double totCharge,double *pot, double dpot[3],double *dmax, double* const Dposition, double* const Dmomentum ){
    /*
    * This function computes the derivative with respect time of position and momenta
    * 
    * Arguments
    * ---------
    * position : vector with the position of the gas particle. *must* be of dimesion = 3
    * momentum : vector with the linear  momentum (that is, velocity * mass) of the gas particle. *must* be of dimesion = 3
    * coord : matrix with the position of each atom. Each row is a different atom an contains the coordinates x,y,z
    * invMass: inverse of the mass. In this context is the inverse of the reduced mass of the system
    * lj : structure with info about the Lennard-Jones parameters
    * charge : vector with the charges of each atom in the molecules
    * dipol : total dipole of the molecule
    * n_atoms : number of atoms within the molecule
    * totCharge : total charge of the molecule
    * 
    * Returns
    * -------
    * pot : total potential energy (LJ + dipole)
    * dpot : vector of dimension 3. derivative of the potential along each direction x,y,z. NOTE: the force is -dpot (*minus* dpot)
    * dmax : minimum distance between the particle of the gas and the molecule (namely, closest atom of the molecule)
    * Dposition : vector with the derivative of the position of the gas particle. Dimension = 3
    * Dmomentum : vector with the derivative of momentum of the gas particle. Dimension = 3
    */
    int i;

    // Derivative of the position = momentum / mass. Literally, we compute momentum * (1/mass). It is faster
    for(i=0;i<3;i++) Dposition[i] = momentum[i] * invMass; 

    // Evaluate the potential
    dljpot(position,coord,lj,charge,dipol,n_atoms,totCharge,pot,dpot,dmax);

    // Derivative of momentum = - Derivative of potential. Note that -derivative of potential is the force.

    for(i=0;i<3;i++) Dmomentum[i] = - dpot[i];

}

void diffeq(double * time,double dt, double* position, double* momentum, double ** coord,const double invMass,const struct LennardJones* const lj,const double* const charge,double dipol, size_t n_atoms,double totCharge,double *pot, double dpot[3],double *dmax, double* const Dposition, double* const Dmomentum  ){
    /*
    * This function integrates the equations of motions using the velocity verlet algorithm.
    * 
    * NOTE: IN/OUT shows the variable that are modified by the function
    * 
    * Arguments
    * ---------
    * time (IN/OUT) : the time of the previous step
    * dt : integration time step
    * position (IN/OUT) : vector with the position of the gas particle. *must* be of dimesion = 3
    * momentum (IN/OUT) : vector with the linear  momentum (that is, velocity * mass) of the gas particle. *must* be of dimesion = 3
    * coord : matrix with the position of each atom. Each row is a different atom an contains the coordinates x,y,z
    * invMass: inverse of the mass. In this context is the inverse of the reduced mass of the system
    * lj : structure with info about the Lennard-Jones parameters
    * charge : vector with the charges of each atom in the molecules
    * dipol : total dipole of the molecule
    * n_atoms : number of atoms within the molecule
    * totCharge : total charge of the molecule
    * Dposition (IN/OUT) : vector with the derivative of the position of the gas particle. *must* be of dimesion = 3
    * Dmomentum (IN/OUT) : vector with the derivative of momentum of the gas particle. *must* be of dimesion = 3
    * 
    * Returns
    * -------
    * time (IN/OUT) : the time of the updated position, velocity and their derivatives
    * pot : total potential energy (LJ + dipole)
    * dpot : vector of dimension 3. derivative of the potential along each direction x,y,z. NOTE: the force is -dpot (*minus* dpot)
    * dmax : minimum distance between the particle of the gas and the molecule (namely, closest atom of the molecule)
    * Dposition (IN/OUT): vector with the derivative of the position of the gas particle. Dimesion = 3
    * Dmomentum (IN/OUT): vector with the derivative of momentum of the gas particle. Dimesion = 3
    * position (IN/OUT) : vector with the position of the gas particle. Dimension = 3
    * momentum (IN/OUT) : vector with the linear  momentum (that is, velocity * mass) of the gas particle. Dimension = 3
    */    
    
    int i;

    // update positions for a full step
    for(i=0;i<3;i++) position[i] += momentum[i]*invMass*dt+0.5*dt*dt*Dmomentum[i]*invMass;

    // update momenta for half step
    for(i=0;i<3;i++) momentum[i] += 0.5*dt*Dmomentum[i];
    
    //update potential and forces
    dljpot(position,coord,lj,charge,dipol,n_atoms,totCharge,pot,dpot,dmax);

    //update the derivative of the momentum
    for(i=0;i<3;i++) Dmomentum[i] = - dpot[i];
    
    // Complete the momentum step    
    for(i=0;i<3;i++) momentum[i] += 0.5*dt*Dmomentum[i];

    // Update the derivative of the position
    for(i=0;i<3;i++) Dposition[i]=momentum[i]*invMass;

    // Update time
    *time=*time+dt;
}

//coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout,ion,parmList):
//return ang,erat,[d1,istep]
void gsang(double ** coord,double theta, double phi, double gamma,const struct LennardJones* const lj,const double* const charge, size_t n_atoms,const struct Molecule *const mol, const double v, const double b, FILE *fp,double * ang, double * erat,int * ifailc){
    /*
    * This function computes the scattering angle using the trajectory method
    * 
    * Arguments
    * ---------
    * coord : coordinates of the atoms within the molecule
    * theta, phi, gamma : angle of the rotation
    * lj : lennard-jones parameters
    * charge : charges of each atom in the molecules
    * n_atoms : number of atoms within the molecule
    * mol : information about the molecule
    * v : velocity of a gas particle
    * b : impact parameter
    * fp : file for logging
    * 
    * Returns
    * -------
    * ang : scattering angle
    * erat : energy ration (kinetic/potential)
    * ifailc : keeps track of how many attempted trajectories are rejected
    * 
    */    

    unsigned int i,ns=0;
    double velocity[3]={0,-v,0},position[3],momentum[3],Dposition[3],Dmomentum[3],top,dt1,dt2,dt,e0,ymin,ymax,pot,dpot[3],dmax,iymin,iymax,iy;
    double etot,time,e,signum;
    
    //printf("pot:%lf dipole:%e charge:%e ep:%e \n",pot,mol->dipole,mol->totalCharge,lj->eolj[0]);

    #if(PRINT_IT) 
    fprintf(fp,"\n specific trajectory parameters\n\n v =% -.4E    b =% -.4E\n",v,b);
    #endif

    // We use a variable time step. So, we should set its value
    top=(v/95.2381)-0.5;
    if (v>=1000) top=10.0;
    if (v>=2000) top=10.0-((v-2000)*7.5e-3);
    if (v>=3000) top=2.5;
    dt1=top*DTSF1*1.0e-11/v;
    dt2=dt1*DTSF2;
    dt=dt1;

    #if(PRINT_IT) 
    fprintf(fp," time steps, dt1 =% -.4E dt2 =% -.4E\n",dt1,dt2);
    #endif

    // We need to locate an appropriate starting position
    position[0]=b;
    position[2]=0.;
    // We start scanning y values between max(y) and min(y)
    e0=0.5*mol->mu*v*v;
    ymax=coord[0][1];
    ymin=coord[0][1];
    for(i=1;i<n_atoms;i++){
        if (ymax<coord[i][1]) ymax=coord[i][1];
        if (ymin>coord[i][1]) ymin=coord[i][1];
    }
    ymax*=1e10;
    ymin*=1e10;

    iymin=trunc(ymin)-1;
    iymax=trunc(ymax)+1;
    iy=iymax;
    position[1]=iy*1e-10;
    //printf("%e %e",iymin,iymax);

    // to make sure that the simulation behaves properly, we should avoid putting the particle
    // too far or too close. We check this, computing the potential at different values of y
    dljpot(position,coord,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
    //printf("pot:%e dipole:%e charge:%e ep:%e y:%e\n",pot,mol->dipole,mol->totalCharge,lj->eolj[0],position[1]);
    if (fabs(pot/e0)<=SW1){
        //printf("iy= %lf first if \n",iy);
        for (iy=iymax-1;iy>iymin;iy--){
            position[1]=iy*1e-10;
            dljpot(position,coord,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
            if (fabs(pot/e0)>=SW1) break;
            if (iy==iymin+1){        
                #if (PRINT_IT) 
                fprintf(fp," trajectory not started - potential too small\n");
                #endif
                *ang=0.0;
                *erat=1.0;
                return;
            }
        }
    }else{
        //printf("iy= %lf second if \n",iy);
        while (fabs(pot/e0)>SW1){
            iy=iy+10;
            position[1]=iy*1e-10;
            dljpot(position,coord,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
        }
        while (fabs(pot/e0)<SW1){
            iy=iy-1;
            position[1]=iy*1e-10;
            dljpot(position,coord,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
        }
    }
    
    // Ok, we are good to go. Let's log a few things first

    // Total energy
    etot=e0+pot;
    
    #if (PRINT_IT) 
    fprintf(fp," trajectory start position =% -.4E\n",iy);
    #endif
    //d1=iy*1e-10;

    // Initial Momentum
    for(i=0;i<3;i++) momentum[i]=velocity[i]*mol->mu;

    // Initial time
    time=0;

    #if (PRINT_IT) 
        fprintf(fp,"\n\n trajectory ns, x,  y,  z,  kin e, dt,    tot e\n");
        fprintf(fp,"                vx, vy, vz, pot e, pot/e0\n\n");
    #endif

    // Initial derivatives
    deriv(position,momentum,coord,mol->invMu,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax,Dposition,Dmomentum);

    // some variables
    *ang=0;
    *erat=0;
    #if (PRINT_IT) 
        // kinetic energy
        e=0.5*mol->invMu*(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]);
        fprintf(fp," %5d % -.4E % -.4E % -.4E % -.4E % -.4E % -.4E\n",ns,position[0],position[1],position[2],e,dt,pot+e);
        fprintf(fp,"       % -.4E % -.4E % -.4E % -.4E % -.4E\n",Dposition[0],Dposition[1],Dposition[2],pot,fabs(pot/e0)); 
    #endif

    // Here starts the molecular dynamics
    while (ns<50){
        do{
            do{
                for(i=0;i<INWR;i++) diffeq(&time,dt,position,momentum,coord,mol->invMu,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax,Dposition,Dmomentum);
                ns+=INWR;
                #if (PRINT_IT)
                    e=0.5*mol->invMu*(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]);
                    fprintf(fp," %5d % -.4E % -.4E % -.4E % -.4E % -.4E % -.4E\n",ns,position[0],position[1],position[2],e,dt,pot+e);
                    fprintf(fp,"       % -.4E % -.4E % -.4E % -.4E % -.4E\n",Dposition[0],Dposition[1],Dposition[2],pot,fabs(pot/e0)); 
                #endif
                if (ns>TRAJ_LOST){
                    fprintf(fp," trajectory lost: b =% -.4E v=% -.4E\n",b,v);
                    *ang=M_PI/2;
                    e=0.5*mol->mu*(Dposition[0]*Dposition[0]+Dposition[1]*Dposition[1]+Dposition[2]*Dposition[2]);
                    *erat=(e+pot)/etot;
                    //istep=ns;
                    return;
                }
                
            }while(dmax<lj->romax);
            if(fabs(pot/e0)>SW2 && dt==dt1) dt=dt2;
            if(fabs(pot/e0)<SW2 && dt==dt2) dt=dt1;
            
        }while(fabs(pot/e0)>SW1);
        
    }

    
    //istep=ns;
    
    // Compute the scattering angle
    if (Dposition[0]>0) signum=1; else if (Dposition[0]<0) signum=-1; else signum=0;
    
    *ang=signum*acos(-Dposition[1]/(sqrt(Dposition[0]*Dposition[0]+Dposition[1]*Dposition[1]+Dposition[2]*Dposition[2])));

    //Check energy conservation
    e=0.5*mol->mu*(Dposition[0]*Dposition[0]+Dposition[1]*Dposition[1]+Dposition[2]*Dposition[2]);
    *erat=(e+pot)/etot;
    if (*erat<1.01 && *erat>0.99) return;  // if everything is ok.otherwise...

    fprintf(fp," \n energy not conserved: e ratio =% -.4E v =% -.4E b =% -.4E\n",*erat,v,b);
    fprintf(fp," gst2 =% -.4E theta =% -.4E phi =% -.4E gamma =% -.4E\n",0.5*mol->mu*v*v/EO,theta*CANG,phi*CANG,gamma*CANG);
    #if(PRINT_IP) 
    fprintf(fp,"\n");
    #endif

    (*ifailc)++;
    if (*ifailc==IFAIL) {
        fclose(fp);
        printf("Energy not conserved");
        exit(EXIT_FAILURE);
    }
}

#if TEST_CTYPES
void rantate(double **InVector,double **OutVector, double * theta, double * phi, double * gamma,size_t n_atoms,double (*RanNum)(void)){
#else
void rantate(double **InVector,double **OutVector, double * theta, double * phi, double * gamma,size_t n_atoms,double (*RanNum)(const gsl_rng *r),const gsl_rng *r){
#endif
    /*
    * This function rotates the molecule under investigationby a random orientation
    * 
    * Arguments
    * ---------
    * InVector : coordinates of the molecules to be rotated    
    * n_atoms : number of atoms within the molecule
    * fp : file for logging data
    * RanNum,g : function for the random number generator. 
    *          Make sure that **the numbers generated are between 0 and 1** (inclusion/exclusion of 1 and/or 0 does not matter). 
    *          As there are many opinions about
    *          the ultimate random number generator, I prefer to use a pointer to function,
    *          instead of hard-coding it. 
    * 
    * 
    * Returns
    * -------
    * OutVector : rotated molecule
    * theta,phi,gamma : (random) angles of the rotation
    */    

    double rnt,rnp,rng;
    
    #if TEST_CTYPES
    rnt=(*RanNum)();
    rnp=(*RanNum)();
    rng=(*RanNum)();
    #else
    rnt=(*RanNum)(r);
    rnp=(*RanNum)(r);
    rng=(*RanNum)(r);
    #endif

    *theta=rnt*2.0*M_PI;
    *phi=asin((rnp*2.0)-1.0)+(M_PI/2.0);
    *gamma=rng*2.0*M_PI;

    Rotate(InVector, OutVector, *theta, *phi, *gamma, n_atoms);

}

void * Thread_MD(void *arg){
    unsigned int ic,ig,im;
    double theta,phi,gamma,**newVector,temp1,temp2,rnb,bst2;
    double b,ang,erat,hold1,hold2,invImp;
    char Filename[20];
    FILE *fp;
    struct arg_MD *args=(struct arg_MD *)arg;

    invImp=1./((double)IMP);
    // Each thread has its own output file
    sprintf(Filename,"thread_%d",args->t_id);
    fp=fopen(Filename,"w");
    if (fp==NULL){
        perror("Error: output file problem ");
        exit(EXIT_FAILURE);
    }

    newVector=AllocDoubleMat(args->n_atoms,3);
    for(ic=0;ic<ITN;ic++){
        
        #if (PRINT_IP) 
        fprintf(fp,"\n cycle number, ic =%3d\n",ic+1);
        #endif

        for(ig=0;ig<INP;ig++){
            #if (PRINT_IP) 
            fprintf(fp,"\n ic =%3d ig =%4d gst2 =% -.4E v =% -.4E\n",ic+1,ig+1,gst2[ig],vList[ig]);
            #endif
            temp1=0.;
            temp2=0.;
            #if (PRINT_IP) 
            fprintf(fp,"\n     b/A        ang      (1-cosX)    e ratio    theta       phi       gamma\n");
            #endif
            for(im=0;im<IMP;im++){
                
                #if TEST_CTYPES
                rnb=(*args->RanNum)();
                rantate(args->coord,newVector, &theta, &phi, &gamma,args->n_atoms,args->RanNum);
                #else
                rnb=(*args->RanNum)(args->r);
                rantate(args->coord,newVector, &theta, &phi, &gamma,args->n_atoms,args->RanNum,args->r);
                #endif

                bst2=rnb*args->b2max[ig];
                b=RO*sqrt(bst2);
                
                #if PRINT_IGS
                    sprintf(Filename,"hold_thread_%d",args->t_id);
                    ftemp=fopen(Filename,"w");
                    fprintf(ftemp,"%d %d %d %d\n",iic,ic+1,ig+1,im+1);
                    fprintf(ftemp,"%lf %lf\n",vList[ig],b);
                    fprintf(ftemp,"%lf %lf %lf\n", theta*CANG,phi*CANG,gamma*CANG);
                    fclose(ftemp);                
                #endif

                gsang(newVector,theta, phi, gamma, args->lj, args->charge, args->n_atoms, args->mol, args->vList[ig], b, fp, &ang, &erat,args->ifailc);
                hold1=1.0-cos(ang);
                hold2=sin(ang)*sin(ang);                
                #if (PRINT_IP) 
                fprintf(fp," % -.4E% -.4E% -.4E% -.4E% -.4E% -.4E% -.4E\n",b*1.e10,ang*CANG,hold1,erat,theta*CANG,phi*CANG,gamma*CANG);
                #endif
                temp1+=(hold1*args->b2max[ig]*invImp);
                temp2+=(1.5*hold2*args->b2max[ig]*invImp);
            }
            args->om11st[ic]+=(temp1*args->wgst[ig]);
            args->om12st[ic]+=(temp1*args->pgst[ig]*args->pgst[ig]*args->wgst[ig]*(1.0/(3.0*args->tst)));
            args->om13st[ic]+=(temp1*(args->pgst[ig]*args->pgst[ig]*args->pgst[ig]*args->pgst[ig])*args->wgst[ig]*(1.0/(12.0*args->tst*args->tst)));
            args->om22st[ic]+=(temp2*args->pgst[ig]*args->pgst[ig]*args->wgst[ig]*(1.0/(3.0*args->tst)));
            args->q1st[ig]+=temp1;
            args->q2st[ig]+=temp2;
            #if (PRINT_IP) 
            fprintf(fp,"\n v =% -.4E     q1st =% -.4E\n\n",vList[ig],q1st[ig]);
            #endif
        }
        #if (PRINT_IP) 
        fprintf(fp,"\n OMEGA(1,1)*=% -.4E\n\n",om11st[ic]);
        #endif
    }
    FreeDoubleMat(newVector);
    fclose(fp);
    return NULL;
}

#if TEST_CTYPES
void mobil2(double **coord,const struct LennardJones* const lj,const double* const charge, size_t n_atoms,const struct Molecule *const mol, unsigned int iic, double *mob,double *cs,double *sdevpc,double (*RanNum)(void),FILE *fp,int print_im2, int *ifailc){
#else
void mobil2(double **coord,const struct LennardJones* const lj,const double* const charge, size_t n_atoms,const struct Molecule *const mol, unsigned int iic, double *mob,double *cs,double *sdevpc,double (*RanNum)(const gsl_rng *r),const gsl_rng *r,FILE *fp,int print_im2, int *ifailc){
#endif
    /*
    * This function computes the cross section and mobility using the trajectory method
    * 
    * Arguments
    * ---------
    * coord : coordinates of the atoms within the molecule
    * lj : lennard-jones parameters
    * charge : charges of each atom in the molecules
    * n_atoms : number of atoms within the molecule
    * mol : information about the molecule
    * iic : an index used to store information about the iteration
    * RanNum,g : function for the random number generator. 
    *          Make sure that **the numbers generated are between 0 and 1** (inclusion/exclusion of 1 and/or 0 does not matter). 
    *          As there are many opinions about
    *          the ultimate random number generator, I prefer to use a pointer to function,
    *          instead of hard-coding it. 
    * 
    * Returns
    * -------
    * mob : mobility
    * cs : cross section
    * sdevpc : standard deviation of OMEGA*(1,1)
    * 
    */    

    unsigned int i;
    int ihold=0,irn,ir,ibst,ig,ic,icc;
    double vector[3]={0,0,0},rotVector[3]={0,0,0},position[3]={0.,0.,0.},rmax=0.,tmp,phi,theta,gamma,rzy,rxy,hold,tst,tst3,gst,dgst;
    double ddd,emaxx,emaxy,emaxz,rmaxx=0,rmaxy=0,rmaxz=0,pot,dmax,dpot[3],r00x=0,r00y=0,r00z=0,sum,sum1,sum2,invSum1;
    double **newVector,wgst[INP],pgst[INP],hold1,hold2,hold3,gstt,invTst,cmin,dbst2,dbst22,gst2[INP],vList[INP],b2max[INP],b,ang,erat,cosx[IBST_MAX],bst2;
    double q1st[N_CORES][INP],q2st[N_CORES][INP],om11st[N_CORES][ITN],om12st[N_CORES][ITN],om13st[N_CORES][ITN],om22st[N_CORES][ITN],mom11st=0.,mom12st=0.,mom13st=0.,mom22st=0.,sdom11st;
    double sterr,ayst,best,cest,term,u2,w,delta,f;
    struct arg_MD args[N_CORES];
    pthread_t tid[N_CORES];
    double t_q1st;

    #if PRINT_IGS
    FILE *ftemp;
    #endif

    newVector=AllocDoubleMat(n_atoms,3);

    if(!print_im2){
        fprintf(fp,"\n mobility calculation by MOBIL2 (trajectory method)\n\n");
        fprintf(fp," global trajectory parameters\n\n sw1 =% -.4E       sw2 =% -.4E\n",SW1,SW2);
        fprintf(fp," dtsf1 =% -.4E     dtsf2 =% -.4E\n",DTSF1,DTSF2);
        fprintf(fp," inwr =%3d              ifail =%5d\n",INWR,IFAIL);
    }

    
    // determine the atom further from the center of mass along x
    if(!print_im2) fprintf(fp,"\n maximum extent orientated along x axis\n");

    for(i=0;i<n_atoms;i++){
        tmp=sqrt(coord[i][0]*coord[i][0]+coord[i][1]*coord[i][1]+coord[i][2]*coord[i][2]);
        //printf("%d %e %e %e %e\n",i,coord[i][0],coord[i][1],coord[i][2],tmp);
        if (tmp>rmax){
            rmax=tmp;
            ihold=i;
        }
    }
    
    for(i=0;i<3;i++) vector[i]=coord[ihold][i]; // store the coordinate of the atom we just found
    
    // Now I use the coordinates of this atom to orient the entire molecules.
    // I position the molecule so that the vector that I just found is oriented along y
    // Notice that the gas particle is "shot" along y

    rzy=sqrt(vector[2]*vector[2]+vector[1]*vector[1]);
    phi=acos(vector[2]/rzy);
    phi+=M_PI_2;
    if (vector[1]<0) phi=(2*M_PI)-phi;
    phi=(2.*M_PI)-phi;
    theta=0.0; 
    gamma=0.0;
    
    Rotate(coord,newVector,theta,phi,gamma, n_atoms);

    #if (PRINT_IU3)
    fprintf(fp,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
    #endif

    for(i=0;i<3;i++) rotVector[i]=newVector[ihold][i]; // keep track of the rotated vector
    rxy=sqrt(rotVector[0]*rotVector[0]+rotVector[1]*rotVector[1]);
    gamma=acos(rotVector[0]/rxy);
    if (rotVector[1]<0) gamma=(2.*M_PI)-gamma;
    gamma=(2.*M_PI)-gamma;

    Rotate(coord,newVector,theta,phi,gamma, n_atoms);

    #if (PRINT_IP==0)
    if (!print_im2) fprintf(fp,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
    #else    
        fprintf(fp,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
        fprintf(fp,"         initial coordinates                        new coordinates\n\n");
        for(i=0;i<n_atoms;i++) fprintf(fp," % -.4E % -.4E % -.4E      % -.4E % -.4E % -.4E\n",coord[i][0],coord[i][1],coord[i][2],newVector[i][0],newVector[i][1],newVector[i][2]);
    #endif
    
    for(i=0;i<3;i++) rotVector[i]=newVector[ihold][i]; // keep track of the rotated vector
    hold=rotVector[0]/rmax;

    if ((hold<0.9999999999) || (hold>1.0000000001) || (rotVector[1]>1.0e-20) || (rotVector[2]>1.0e-20) || (rotVector[1]<-1.0e-20) || (rotVector[2]<-1.0e-20)){
        fprintf(fp,"\n Problem orientating along x axis\n\n");
        for(i=0;i<n_atoms;i++) fprintf(fp,"% -.4E % -.4E % -.4E % -.4E\n",newVector[i][0],newVector[i][1],newVector[i][2],sqrt(newVector[i][0]*newVector[i][0]+newVector[i][1]*newVector[i][1]+newVector[i][2]*newVector[i][2]));
    }

    #if (PRINT_IP) 
    fprintf(fp,"\n\n");
    #endif

    irn=1000;  // TODO : Add to constant list?
    ddd=(rmax+lj->romax)/(double)irn;

    for(i=0;i<3;i++) position[i]=0;
    emaxx=0.0;
    for(ir=0;ir<irn;ir++){
        position[0]=rmax+lj->romax-((double)ir*ddd);
        dljpot(position,newVector,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
        if (pot>0) break;
        r00x=position[0];
        if  (pot<emaxx){
            rmaxx=position[0];
            emaxx=pot;
        }
    }
    if (!print_im2) fprintf(fp," along x axis emax =% -.4EeV rmax =% -.4EA r00 =% -.4EA\n",emaxx/XE,rmaxx*1.0e10,r00x*1.0e10);

    for(i=0;i<3;i++) position[i]=0;
    emaxy=0.0;
    for(ir=0;ir<irn;ir++){
        position[1]=rmax+lj->romax-((double)ir*ddd);
        dljpot(position,newVector,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
        if (pot>0) break;
        r00y=position[1];
        if  (pot<emaxy){
            rmaxy=position[1];
            emaxy=pot;
        }
    }
    if (!print_im2) fprintf(fp," along y axis emax =% -.4EeV rmax =% -.4EA r00 =% -.4EA\n",emaxy/XE,rmaxy*1.0e10,r00y*1.0e10);

    for(i=0;i<3;i++) position[i]=0;
    emaxz=0.0;
    for(ir=0;ir<irn;ir++){
        position[2]=rmax+lj->romax-((double)ir*ddd);
        dljpot(position,newVector,lj,charge,mol->dipole,n_atoms,mol->totalCharge,&pot,dpot,&dmax);
        if (pot>0) break;
        r00z=position[2];
        if  (pot<emaxz){
            rmaxz=position[2];
            emaxz=pot;
        }
    }
    if (!print_im2) fprintf(fp," along z axis emax =% -.4EeV rmax =% -.4EA r00 =% -.4EA\n\n",emaxz/XE,rmaxz*1.0e10,r00z*1.0e10);

    // Set up integration over gst
    tst=XK*mol->temperature/EO;
    invTst=1/tst;
    if (!print_im2) fprintf(fp,"\n t*=% -.4E\n",tst);
    tst3=tst*tst*tst;
    dgst=5e-7*6.*sqrt(tst);
    gst=dgst;
    sum=0.0;
    sum1=0.0;
    sum2=0.0;
    
    if (!print_im2){
        fprintf(fp,"\n\n set-up gst integration - integration over velocity\n");
        fprintf(fp,"\n     pgst        wgst         v         ke/kt       gst^5*     frac of\n");
        fprintf(fp,"                                                exp(gst^2/tst)   sum\n\n");
    }

    for(i=1;i<INP+1;i++) sum1 +=sqrt((double)i);
    invSum1=1./sum1;

    for(i=0;i<INP;i++){
        hold1=sqrt((double)i+1.);
        wgst[i]=hold1*invSum1;
        sum2+=sqrt((double)i);
        gstt=tst3*(sum2+(0.5*hold1))*invSum1;

        while(1){   
            sum+=(exp(-gst*gst*invTst)*gst*gst*gst*gst*gst*dgst);
            gst=gst+dgst;
            if (sum>gstt){
                pgst[i]=gst-(0.5*dgst);
                break;
            }
            if (sum==gstt) break;
        }

        if (!print_im2){
            hold1=sqrt((pgst[i]*pgst[i]*EO)/(0.5*mol->mu));
            hold2=0.5*mol->mu*hold1*hold1/(XK*mol->temperature);
            hold3=exp(-pgst[i]*pgst[i]*invTst)*pgst[i]*pgst[i]*pgst[i]*pgst[i]*pgst[i];
            fprintf(fp," % -.4E % -.4E % -.4E % -.4E % -.4E % -.4E\n",pgst[i],wgst[i],hold1,hold2,hold3,sum/tst3);
        }
        gst2[i]=pgst[i]*pgst[i];
        vList[i]=sqrt((gst2[i]*EO)/(0.5*mol->mu));    
    }
    // Let's find b2max next
    
    dbst2=1.0;
    dbst22=0.1*dbst2;
    cmin=0.0005;
    if (!print_im2) fprintf(fp,"\n\n set up b2 integration - integration over impact parameter\n\n minimum value of (1-cosX) =% -.4E\n\n",cmin);

    for (ig=INP-1;ig>-1;ig--){
        ibst= (int) trunc(rmaxx/RO)-6;
        if (ig<INP-1) ibst=(int) trunc(b2max[ig+1]/dbst2)-6;
        if (ibst<0) ibst=0;

        #if(PRINT_IP)
            fprintf(fp,"\n gst2 =% -.4E v =% -.4E\n",gst2[ig],vList[ig]);
            fprintf(fp,"      b          bst2       X ang       cos(X)      e ratio\n");
        #endif
        
        while (1){        
            bst2=dbst2*(double) ibst;
            b=RO*sqrt(bst2);
            gsang(newVector,theta, phi, gamma, lj, charge, n_atoms, mol, vList[ig], b, fp,&ang, &erat,ifailc);
            cosx[ibst]=1.0-cos(ang);

            #if(PRINT_IP) 
            fprintf(fp," % -.4E % -.4E % -.4E % -.4E % -.4E\n",b,bst2,ang,cosx[ibst],erat);
            #endif

            if (ibst<5  || (ibst>=5 && (!(cosx[ibst]<cmin && cosx[ibst-1]<cmin && cosx[ibst-2] < cmin && cosx[ibst-3]<cmin && cosx[ibst-4]<cmin)))){
                ibst++;
                if (ibst>=IBST_MAX){
                    printf(" ibst greater than %d",IBST_MAX);
                    fclose(fp);
                    exit(EXIT_FAILURE);
                }
            }else break;
        }

        b2max[ig]=(double)(ibst-5)*dbst2;
        do {
            b2max[ig]=b2max[ig]+dbst22;
            b=RO*sqrt(b2max[ig]);
            gsang(newVector,theta, phi, gamma, lj, charge, n_atoms, mol, vList[ig], b, fp,&ang, &erat,ifailc);
        }while(1.0-cos(ang)>cmin);
    }   
    if(!print_im2){
        fprintf(fp,"\n     gst           b2max/ro2         b/A\n\n");
        for (ig=0;ig<INP;ig++)  fprintf(fp," % -.4E     % -.4E     % -.4E\n",pgst[ig],b2max[ig],RO*sqrt(b2max[ig])*1.0e10);
    }
    
    //  Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
    //  by integrating Q(1)* or Q(2)* over all orientations, and initial 
    //  relative velocities.


    if(!print_im2){
        fprintf(fp,"\n\n number of complete cycles (itn) =%6d\n",ITN*N_CORES);
        fprintf(fp," number of velocity points (inp) =%6d\n",INP);
        fprintf(fp," number of random points (imp) =%6d\n",IMP);
        fprintf(fp," total number of points =%7d\n\n",ITN*INP*IMP*N_CORES);
    }
    #if (PRINT_IP) 
    fprintf(fp," start mobility calculation\n");
    #endif
    
    memset(q1st,0,INP*N_CORES*sizeof(double));
    memset(q2st,0,INP*N_CORES*sizeof(double));
    memset(om11st,0,ITN*N_CORES*sizeof(double));
    memset(om12st,0,ITN*N_CORES*sizeof(double));
    memset(om13st,0,ITN*N_CORES*sizeof(double));
    memset(om22st,0,ITN*N_CORES*sizeof(double));

    //printf("OK-a\n");
    for(i=0;i<N_CORES;i++){
        SetArgMD(&args[i],i,coord,lj,charge,mol,RanNum,r, n_atoms,b2max, vList,wgst,pgst,ifailc,tst,om11st[i],om12st[i],om13st[i],om22st[i],q1st[i],q2st[i]);
        pthread_create(&tid[i],NULL,Thread_MD,(void *) &args[i]);
        //pthread_join(tid[i],NULL);
    }
    for(i=0;i<N_CORES;i++) pthread_join(tid[i],NULL);
    //printf("OK-b\n");

    // for(ic=0;ic<ITN;ic++){
        
    //     #if (PRINT_IP) 
    //     fprintf(fp,"\n cycle number, ic =%3d\n",ic+1);
    //     #endif

    //     for(ig=0;ig<INP;ig++){
    //         #if (PRINT_IP) 
    //         fprintf(fp,"\n ic =%3d ig =%4d gst2 =% -.4E v =% -.4E\n",ic+1,ig+1,gst2[ig],vList[ig]);
    //         #endif
    //         temp1=0.;
    //         temp2=0.;
    //         #if (PRINT_IP) 
    //         fprintf(fp,"\n     b/A        ang      (1-cosX)    e ratio    theta       phi       gamma\n");
    //         #endif
    //         for(im=0;im<IMP;im++){
                
    //             #if TEST_CTYPES
    //             rnb=(*RanNum)();
    //             rantate(coord,newVector, &theta, &phi, &gamma,n_atoms,RanNum);
    //             #else
    //             rnb=(*RanNum)(r);
    //             rantate(coord,newVector, &theta, &phi, &gamma,n_atoms,RanNum,r);
    //             #endif

    //             bst2=rnb*b2max[ig];
    //             b=RO*sqrt(bst2);
                
    //             #if PRINT_IGS
    //                 ftemp=fopen("hold_c","w");
    //                 fprintf(ftemp,"%d %d %d %d\n",iic,ic+1,ig+1,im+1);
    //                 fprintf(ftemp,"%lf %lf\n",vList[ig],b);
    //                 fprintf(ftemp,"%lf %lf %lf\n", theta*CANG,phi*CANG,gamma*CANG);
    //                 fclose(ftemp);
                
    //             #endif

    //             gsang(newVector,theta, phi, gamma, lj, charge, n_atoms, mol, vList[ig], b, fp, &ang, &erat,ifailc);
    //             hold1=1.0-cos(ang);
    //             hold2=sin(ang)*sin(ang);                
    //             #if (PRINT_IP) 
    //             fprintf(fp," % -.4E% -.4E% -.4E% -.4E% -.4E% -.4E% -.4E\n",b*1.e10,ang*CANG,hold1,erat,theta*CANG,phi*CANG,gamma*CANG);
    //             #endif
    //             temp1+=(hold1*b2max[ig]*invImp);
    //             temp2+=(1.5*hold2*b2max[ig]*invImp);
    //         }
    //         om11st[ic]+=(temp1*wgst[ig]);
    //         om12st[ic]+=(temp1*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)));
    //         om13st[ic]+=(temp1*(pgst[ig]*pgst[ig]*pgst[ig]*pgst[ig])*wgst[ig]*(1.0/(12.0*tst*tst)));
    //         om22st[ic]+=(temp2*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)));
    //         q1st[ig]+=temp1;
    //         q2st[ig]+=temp2;
    //         #if (PRINT_IP) 
    //         fprintf(fp,"\n v =% -.4E     q1st =% -.4E\n\n",vList[ig],q1st[ig]);
    //         #endif
    //     }
    //     #if (PRINT_IP) 
    //     fprintf(fp,"\n OMEGA(1,1)*=% -.4E\n\n",om11st[ic]);
    //     #endif
    // }

    //     calculate running averages

    hold1=0.;
    hold2=0.;
    if(!print_im2) fprintf(fp,"\n summary of mobility calculations\n\n cycle     cs/A^2      avge cs/A^2        Ko^-1       avge Ko^-1\n");
    
    for(i=0;i<N_CORES;i++){
        for(icc=0;icc<ITN;icc++){
            tmp=1.0/(mol->mobility/(sqrt(mol->temperature)*om11st[i][icc]*M_PI*RO2));
            hold1+=om11st[i][icc];
            hold2+=tmp;
            if(!print_im2) fprintf(fp," %3d    % -.4E    % -.4E    % -.4E    % -.4E\n",(i*ITN)+icc+1,om11st[i][icc]*M_PI*RO2*1.e20,hold1*M_PI*RO2*1e20/((double)((i*ITN)+icc+1)),tmp,hold2/((double)((i*ITN)+icc+1)));        
        }
    }
    if(!print_im2) fprintf(fp,"\n\n average values for q1st\n\n     gst2        wgst        q1st\n");

    if(!print_im2) for(ig=0;ig<INP;ig++) {
        t_q1st=0.;
        for(i=0;i<N_CORES;i++) t_q1st+=q1st[i][ig];
        fprintf(fp," % -.4E % -.4E % -.4E\n",pgst[ig]*pgst[ig],wgst[ig],t_q1st/(double) (INP));
    }
    for(i=0;i<N_CORES;i++) for(ic=0;ic<ITN;ic++) mom11st+=om11st[i][ic];
    for(i=0;i<N_CORES;i++) for(ic=0;ic<ITN;ic++) mom12st+=om12st[i][ic];
    for(i=0;i<N_CORES;i++) for(ic=0;ic<ITN;ic++) mom13st+=om13st[i][ic];
    for(i=0;i<N_CORES;i++) for(ic=0;ic<ITN;ic++) mom22st+=om22st[i][ic];
    mom11st/=(double) (ITN*N_CORES);
    mom12st/=(double) (ITN*N_CORES);
    mom13st/=(double) (ITN*N_CORES);
    mom22st/=(double) (ITN*N_CORES);
    //standard deviation
    sdom11st=0.0;
    for(i=0;i<N_CORES;i++) for(ic=0;ic<ITN;ic++) sdom11st+=(om11st[i][ic]-mom11st)*(om11st[i][ic]-mom11st);
    sdom11st=sqrt(sdom11st/(double) (ITN*N_CORES));
    sterr=sdom11st/sqrt((double)(ITN*N_CORES));

    if(!print_im2) fprintf(fp,"\n\n mean OMEGA*(1,1) =% -.4E\n standard deviation =% -.4E\n standard error of mean =% -.4E\n",mom11st,sdom11st,sterr);
    *cs=mom11st*M_PI*RO2;
    *sdevpc=100.0*sdom11st/mom11st;

    //Use omegas to obtain higher order correction factor to mobility

    ayst=mom22st/mom11st;
    best=((5.0*mom12st)-(4.0*mom13st))/mom11st;
    cest=mom12st/mom11st;
    term=((4.0*ayst)/(15.0))+(.50*((mol->TotalMass-mol->mGas)*(mol->TotalMass-mol->mGas))/(mol->mGas*mol->TotalMass));
    u2=term-(.083330*(2.40*best+1.0)*(mol->mGas/mol->TotalMass));
    w=(mol->mGas/mol->TotalMass);
    delta=((((6.0*cest)-5.0)*((6.0*cest)-5.0))*w)/(60.0*(1.0+u2));
    f=1.0/(1.0-delta);
    if (!print_im2){
        fprintf(fp,"\n\n f value for second order correction=% -.4E\n",f);
        fprintf(fp," (integrations for second order correction are not\n");
        fprintf(fp," accurate, check if correction becomes significant)\n");
        fprintf(fp,"\n omega*12 =% -.4E  omega*13 =% -.4E  omega*22 =% -.4E\n",mom12st,mom13st,mom22st);
        fprintf(fp,"       u2 =% -.4E         w =% -.4E     delta =% -.4E\n",u2,w,delta);
    }
    *mob=(mol->mobility*f)/(sqrt(mol->temperature)*(*cs));
    fprintf(fp,"\n\n average (second order) TM mobility =% -.4E\n",*mob);
    fprintf(fp," inverse average (second order) TM mobility =% -.4E\n",1.0/(*mob));
    fprintf(fp," average TM cross section =% -.4E\n",(*cs)*1.e20);

    //free a bounch of stuff
    FreeDoubleMat(newVector);
}


void che(int im,double **coord,struct HardSphere* hs,double *cof,double cop,double *versor, double *yr, double *zr, int *kp,size_t n_atoms){
    /*
    * This function checks whether the gas particle hits the molecule or is flying away
    * 
    * Arguments
    * ---------
    * im : order of reflection/collision (== how many times it bounces against the molecule before flying away)
    * coord : coordinates of the atoms within the molecule
    * hs : hard sphere parameters
    * cop : previous collision cross section
    * yr,zr (IN/OUT) : coordinates of the collision point (IN : used only for the first collision, then both =0)
    * kp (IN/OUT) : there was *already* a collision (1) or not (0). 
    * n_atoms : number of atoms within the molecule
    * 
    * Returns
    * -------
    * cof : collision cross section after a collision
    * yr,zr (IN/OUT) : coordinates of the collision point
    * kp (IN/OUT) : there was a collision (1) or not (0) 
    * 
    */    

    double xl=XL_INIT;  // just a number larger than any cluster
    unsigned int i;
    double ras,dev,xc,xv,yv,zv,xve1,xve2,xve3,yve1,yve2,yve3,zve1,zve2,zve3,xxv,xyv,xzv,xyz,rad2,adr1,adr2,ad1,yd,zd;
    double xne,yne;
    int ki;

    /*
    * This function "follows" the trajectory of the particle *after* the first collision
    * The location of the first collision is computed in the PA part of the algorithm and passed
    * to this function. Orginal code by Alexandre Shvartsburg.
    * 
    * 
    */

    // The PA algorithms aligns the direction *after* the collision with the x-axis.
    // Thus, we do not need to align at first
    // DIFF: im starts from 1 in original code
    if (im!=0){
        *yr=0.;
        *zr=0.;
    }

    ki=-1;
    for(i=0;i<n_atoms;i++){
        if(coord[i][0]>1e-16 || im==0){ // what I mean is that if x is not zero. Being a float x>0 within numeric precision
            // now we check for collisions. Remember that the gas particle arrives along the x axis (or better, the molecule is oriented in that way)
            // yd,zd = collision coordinates, dev=distance between atom and collision point
            yd=*yr-coord[i][1];
            zd=*zr-coord[i][2];
            ras=(yd*yd)+(zd*zd);
            dev=sqrt(ras);

            // Let's check if the collision with atom i actually took place
            // We look for the atom closest (along x) to the collision point

            if (dev<=hs->rhs[i]){
                xc=coord[i][0]-sqrt(hs->rhs2[i]-ras);
                if (xc<xl){
                    xl=xc;
                    ki=i;
                }
            }
        }   

    }
    // If there was a collision, let's compute the distance from the closest atom
    if(ki>=0){
        *kp=1;
        xv=xl-coord[ki][0];
        yv=*yr-coord[ki][1];
        zv=*zr-coord[ki][2];
    

        // Now we perform a translation so that the collison point (xl,yr,zr) becomes (0,0,0)
        for(i=0;i<n_atoms;i++){
            coord[i][0]-=xl;
            coord[i][1]-=*yr;
            coord[i][2]-=*zr;
        }

        // Then, we rotate the molecule so that the *reflected* (namely, after impact) particle
        // is along x. Note that y and z are arbitrary as there is a infinite number of such transformations
        // xve1,xve2,xve3 is the direction vector of the colliding particle
        
        // Lets define the rotation matrix elements
        xxv=2.0*xv*xv;
        xyv=2.0*xv*yv;
        xzv=2.0*xv*zv;
        xyz=xyv*xzv;
        rad2=hs->rhs2[ki]-xxv;
        ad1=(rad2*rad2)+(xyv*xyv);
        adr1=sqrt(ad1);
        adr2=sqrt(ad1*ad1+xyz*xyz+xzv*xzv*rad2*rad2);
        xve1=1.0-2.0*xv*xv/hs->rhs2[ki];
        xve2=-xyv/hs->rhs2[ki];
        xve3=-xzv/hs->rhs2[ki];
        yve1=xyv/adr1;
        yve2=rad2/adr1;
        yve3=0.0;
        zve1=rad2*xzv/adr2;
        zve2=-xyz/adr2;
        zve3=ad1/adr2;

        // Now we rotate
        for(i=0;i<n_atoms;i++){
            xne=xve1*coord[i][0]+xve2*coord[i][1]+xve3*coord[i][2];
            yne=yve1*coord[i][0]+yve2*coord[i][1]+yve3*coord[i][2];
            coord[i][2]=zve1*coord[i][0]+zve2*coord[i][1]+zve3*coord[i][2];
            coord[i][0]=xne;
            coord[i][1]=yne;
        }        
        xne=xve1*versor[0]+xve2*versor[1]+xve3*versor[2];
        yne=yve1*versor[0]+yve2*versor[1]+yve3*versor[2];
        versor[2]=zve1*versor[0]+zve2*versor[1]+zve3*versor[2];
        versor[0]=xne;
        versor[1]=yne;
        *cof=cos(0.50*(M_PI-acos(versor[0])));
    }else *cof=cop;
}


void *thread_MC(void *arg){

    double versor[3],cof[INOR+1],**newVector;
    double theta=0.,phi=0.,gamma=0.,ymin=0.,ymax=0.,zmin=0.,zmax=0.,ydi,zdi,pls,yr,zr;
    int kp=0;
    struct arg_mc *args=(struct arg_mc*)arg;
    unsigned int i,n,im,imn;

    newVector=AllocDoubleMat(args->n_atoms,3);
    *args->imm=1;
    *args->crb=0.;

    for(i=0;i<INUM;i++){
        #if TEST_CTYPES
        rantate(args->coord,newVector, &theta, &phi, &gamma,args->n_atoms,args->RanNum);
        #else
        rantate(args->coord,newVector, &theta, &phi, &gamma,args->n_atoms,args->RanNum,args->r);
        #endif


        #if (PRINT_IU2 || PRINT_IU3)
        fprintf(fp,"\n\n coordinates rotated by ROTATE\n\n theta=% -.4E phi=% -.4E gamma=% -.4E\n\n",theta*CANG,phi*CANG,gamma*CANG);
        #endif
        
        #if (PRINT_IU2)
        fprintf(fp,"         initial coordinates                        new coordinates\n\n");
        for(n=0;n<args->n_atoms;n++) fprintf(fp," % -.4E % -.4E % -.4E      % -.4E % -.4E % -.4E\n",args->coord[n][0],args->coord[n][1],args->coord[n][2],newVector[n][0],newVector[n][1],newVector[n][2]);
        #endif

        // Build  a  box around the molecule (the particle comes from the x direction)
        ymin=0.;ymax=0.;zmin=0.;zmax=0.;
        for(n=0;n<args->n_atoms;n++){
            if (newVector[n][1]-args->hs->rhs[n]<ymin) ymin=newVector[n][1]-args->hs->rhs[n];
            if (newVector[n][1]+args->hs->rhs[n]>ymax) ymax=newVector[n][1]+args->hs->rhs[n];
            if (newVector[n][2]-args->hs->rhs[n]<zmin) zmin=newVector[n][2]-args->hs->rhs[n];
            if (newVector[n][2]+args->hs->rhs[n]>zmax) zmax=newVector[n][2]+args->hs->rhs[n];
        }
        // Compute Area of the box:
        ydi=ymax-ymin;
        zdi=zmax-zmin;
        pls=ydi*zdi;

        //Pick a random point within the box:
        #if TEST_CTYPES
        yr=ymin+ydi*(*args->RanNum)();
        zr=zmin+zdi*(*args->RanNum)();
        #else
        yr=ymin+ydi*(*args->RanNum)(args->r);
        zr=zmin+zdi*(*args->RanNum)(args->r);
        #endif

        cof[0]=0.;
        kp=0;
        versor[0]=1.; versor[1]=0.; versor[2]=0.;
        for(im=0;im<INOR+1;im++){
            //printf("yr:%le zr:%le\n",yr,zr);
            che(im,newVector,args->hs,&cof[im+1],cof[im],versor, &yr, &zr, &kp,args->n_atoms);
            // If there is no collision (kp=0), then the angle before and after "no collision" does not change
            if( cof[im+1] == cof[im] ){
                //printf("%d %lf %lf\n",im,cof[im+1],cof[im]); 
                if(im > *args->imm ) *args->imm = im;
                for (imn=im+2;imn<INOR+1;imn++) cof[imn]=cof[im+1];  // we set all the "collision" to the last actual collision
                break;                
            } else if (im+1 > *args->imm) *args->imm = im+1;  // if we got a collision (kp=1)
        }
        if(kp==1) *args->crb += pls; // Projection approximation (= simple monte carlo = count only first collision)
        for(im=0;im<INOR;im++) args->crof[im] += pls*cof[im+1]*cof[im+1];  // for hard sphere (= count multiple collisions)
    }
    FreeDoubleMat(newVector);
    return NULL;
}


#if TEST_CTYPES
void mobil4(double **coord,struct HardSphere* hs,const struct Molecule *const mol,double (*RanNum)(void),double *ehscs,double *ehsmob,double *pacs,double *pamob,unsigned int * imm,size_t n_atoms,FILE *fp,int print_im4){
#else
void mobil4(double **coord,struct HardSphere* hs,const struct Molecule *const mol,double (*RanNum)(const gsl_rng *r),const gsl_rng *r,double *ehscs,double *ehsmob,double *pacs,double *pamob,unsigned int * imm,size_t n_atoms,FILE *fp,int print_im4){
#endif
    /*
    * This function computes the cross section and mobility using the projection approximation and hard sphere models
    * 
    * Arguments
    * ---------
    * coord : coordinates of the atoms within the molecule
    * hs : hard sphere parameters
    * n_atoms : number of atoms within the molecule
    * mol : information about the molecule
    * print_im4 : printing verbosity ( 0 : maximum verbosity; 1 : minimum verbosity)
    * RanNum,g : function for the random number generator. 
    *          Make sure that **the numbers generated are between 0 and 1** (inclusion/exclusion of 1 and/or 0 does not matter). 
    *          As there are many opinions about
    *          the ultimate random number generator, I prefer to use a pointer to function,
    *          instead of hard-coding it. 
    * fp : file for logging
    * 
    * Returns
    * -------
    * ehsmob,ehscs : mobility and cross section from hard sphere multiple scattering
    * pamob,pacs : mobility and cross section from projection approximation
    * imm : how many collision before leaving the molecule
    * 
    */    

    double crof_t[N_CORES][INOR],crb=0.,crof[INOR];
    double crb_t[N_CORES];
    unsigned int i,imm_t[N_CORES],im;
    pthread_t tid[N_CORES];
    struct arg_mc args[N_CORES];

    memset(crof,0,sizeof(double)*INOR);
    memset(crof_t,0,sizeof(double)*INOR*N_CORES);

    *imm=1;

    // BEGIN MONTE CARLO
    for(i=0;i<N_CORES;i++){
        #if TEST_CTYPES
        SetArgmc(&args[i],coord,hs,RanNum,&crb_t[i],&imm_t[i],crof_t[i],n_atoms);
        #else
        SetArgmc(&args[i],coord,hs,RanNum,r,&crb_t[i],&imm_t[i],crof_t[i],n_atoms);
        #endif
        pthread_create(&tid[i],NULL,&thread_MC,(void *) &args[i]);
    }
    for(i=0;i<N_CORES;i++) pthread_join(tid[i],NULL);

    // END of MONTE CARLO
    for(i=0;i<N_CORES;i++) if(imm_t[i]>*imm) *imm=imm_t[i];
    for(i=0;i<N_CORES;i++) crb+=crb_t[i];
    for(im=0;im<INOR;im++) for(i=0;i<N_CORES;i++) crof[im]+=crof_t[i][im];
    
    *pacs=crb/((double)(INUM*N_CORES));
    *pamob=mol->mobility/(sqrt(mol->temperature)*(*pacs));

    for(im=0;im<INOR;im++) crof[im] /=((double)(INUM*N_CORES)*0.5);
    *ehscs=crof[*imm-1];
    *ehsmob=mol->mobility/(sqrt(mol->temperature)*(*ehscs));

    // Output results
    if (print_im4==0){
        fprintf(fp,"\n\n mobility calculation by MOBIL4 (HS scattering)\n");
        fprintf(fp,"\n number of Monte Carlo trajectories =%7d\n",INUM*N_CORES);
        fprintf(fp," maximum number of reflections allowed =%3d\n",INOR);
    }

    fprintf(fp,"\n average PA mobility =% -.4E\n",*pamob);
    fprintf(fp," inverse average PA mobility =% -.4E\n",1.0/(*pamob));
    fprintf(fp," average PA cross section =% -.4E\n",(*pacs)*1.e20);
    fprintf(fp,"\n maximum number of reflections encountered =%3d\n",*imm);

    if (print_im4==0){
        fprintf(fp,"\n     order   cross section\n");
        for (im=1;im<=*imm;im++) fprintf(fp,"     %3d     % -.5E\n",im,crof[im-1]);
    }

    fprintf(fp,"\n average EHS mobility =% -.4E\n",*ehsmob);
    fprintf(fp," inverse average EHS mobility =% -.4E\n",1.0/(*ehsmob));
    fprintf(fp," average EHS cross section =% -.4E\n",(*ehscs)*1.e20);
}

#if TEST_CTYPES
void wrap_mobil2(double **coord,const struct LennardJones* const lj,const double* const charge, size_t n_atoms,const struct Molecule *const mol, unsigned int iic, double *mob,double *cs,double *sdevpc,double (*RanNum)(void),int print_im2){

    FILE * fp;
    int ifailc=0;
    fp=fopen("out.c.dat","w");

    mobil2(coord, lj,charge, n_atoms,mol, iic, mob, cs, sdevpc, *RanNum, fp, print_im2, &ifailc);
    fclose(fp);
}

void wrap_mobil4(double **coord,struct HardSphere* hs,const struct Molecule *const mol,double (*RanNum)(void),double *ehscs,double *ehsmob,double *pacs,double *pamob,unsigned int * imm,size_t n_atoms){
    /*
    * mobil4 wrapper to be used with ctypes.
    * ctypes does not share files, so the wrapper open/closes it independently
    */
   
    FILE *fp;
    fp=fopen("out.c.dat","w");
    mobil4(coord, hs, mol,RanNum,ehscs,ehsmob,pacs,pamob,imm,n_atoms,fp,0);
    fclose(fp);
}
#endif

#if TEST_CTYPES
int wrap_main(double (*RanNum)(void)){
#else
int main(void){
    double (*RanNum)(const gsl_rng *);
    gsl_rng * r;
    unsigned long int Rseed=0;
#endif
    struct LennardJones lj;
    struct Molecule mol;
    struct inputFileData inData;    
    struct AtomData atoms;
    struct HardSphere hs;
    
    double asymp,*tmm,*tmc,*asympp,*ehsc,*ehsm,*pac,*pam,sdevpc=0.;
    double pacs=0.,pamob=0.,ehscs=0.,ehsmob=0.,aasymp=0.,tmcs=0.,tmmob=0.;
    unsigned int imm,immmin,immmax,iic,i;
    int im2=0,im4=0,ifailc=0;
    FILE *fpIn,*fpOut;

    #if TEST_CTYPES
    char * inFile="a10A1.mfj";
    char * outFile="out.c.data";
    
    #else
    char buffer[MAX_LINE],t_buffer[MAX_LINE],inFile[MAX_LINE],outFile[MAX_LINE];
    char * runFile="mobcal.run";    
    fpIn=fopen(runFile,"r");
    if (fpIn==NULL){
        perror("Error: Run file problem ");
        exit(EXIT_FAILURE); 
    }

    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        perror("Error: Input file name was not found ");
        exit(EXIT_FAILURE);
    }else {
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
        strcpy(t_buffer,strtok(buffer," \t\r\n"));  // remove trailing non-characters
        sscanf(t_buffer,"%s",inFile);
    }
    
    if (fgets(buffer,MAX_LINE,fpIn)==NULL){
        fclose(fpIn);
        perror("Error: Output file name was not found ");
        exit(EXIT_FAILURE);
    }else {
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
        strcpy(t_buffer,strtok(buffer," \t\r\n"));  // remove trailing non-characters
        sscanf(t_buffer,"%s",outFile);
    }

    if (fgets(buffer,MAX_LINE,fpIn)!=NULL){
        if (buffer[strlen(buffer)-1]!='\n') while(getc(fpIn)!='\n'){}; //consume rest of the line (if any)
        sscanf(buffer,"%lu",&Rseed);
    }else Rseed=0; // NOTE: in gsl =0 means "use default seed". This seed is different for different algorithms and can be changed by environment. Check docs

    fclose(fpIn);
    #endif

    fpIn=fopen(inFile,"r");
    if (fpIn==NULL){
        perror("Error: input file with structure problem ");
        exit(EXIT_FAILURE);
    }
    
    fpOut=fopen(outFile,"w");
    if (fpOut==NULL){
        fclose(fpIn);
        perror("Error: output file problem ");
        exit(EXIT_FAILURE);
    }

    // Read the input file and set up all the data structures
    fcoord(inFile,&inData,&lj,&hs,&mol,&atoms,&asymp,fpIn,fpOut);

    // Print some constants
    fprintf(fpOut," Lennard-Jones scaling parameters: eo=% -.4E ro=% -.4E\n",EO/XE,RO*1e10);

    // Alloc some arrays
    tmm=AllocDoubleVec(inData.n_coord);
    tmc=AllocDoubleVec(inData.n_coord);
    asympp=AllocDoubleVec(inData.n_coord);
    ehsc=AllocDoubleVec(inData.n_coord);
    ehsm=AllocDoubleVec(inData.n_coord);
    pac=AllocDoubleVec(inData.n_coord);
    pam=AllocDoubleVec(inData.n_coord);
    asympp[0]=asymp;
    
    //  ion-induced dipole potential with Helium (supposed to be the buffer gas)
    //  Polarizability of helium = 0.204956d-30 m^3
    //  xeo is permitivity of vacuum, 8.854187817d-12 Farad/meter 

    mol.dipole=XE*XE*0.204956e-30/(8*M_PI*XEO);
    fprintf(fpOut," dipole constant =% -.4E\n",mol.dipole);

    // Set the temperature
    mol.temperature=TEMPERATURE; 

    // Set the mass of the buffer (Helium)
    mol.mGas=4.0026;
    mol.mu=((mol.TotalMass*mol.mGas)/(mol.TotalMass+mol.mGas))/(XN*1000);
    mol.invMu=1./mol.mu;

    // Here we evaluate the mobility constant
    mol.mobility=sqrt(18*M_PI)/16;
    mol.mobility*=sqrt(XN*1000)*sqrt((1./mol.mGas)+(1./mol.TotalMass));
    mol.mobility*=XE/sqrt(XK);
    mol.mobility /=(XN/XMV);
    fprintf(fpOut," mobility constant =% -.4E\n",mol.mobility); 
    fprintf(fpOut," temperature =% -.4E\n",mol.temperature);

    // DIFF: In the original code follows an initialization of the random number generator.
    // Here I am using a library  defined algorithm
    #if !TEST_CTYPES
    r=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,Rseed);
    RanNum=&gsl_rng_uniform; 
    #endif

    // Some bookkeping for later use.
    imm=0;
    immmin=INOR;
    immmax=0;

#if TEST
    inData.n_coord=2;
#endif

    for (iic=0;iic<inData.n_coord;iic++){
        fprintf(fpOut,"\n\n");
        if (inData.n_coord > 1) fprintf(fpOut,"\n coordinate set =%5d\n",iic+1);
        fprintf(fpOut,"\n structural asymmetry parameter =%8.4f\n",asympp[iic]);

        #if TEST_CTYPES
        mobil4(atoms.coord, &hs, &mol,RanNum,&ehsc[iic],&ehsm[iic],&pac[iic],&pam[iic],&imm,inData.n_atoms,fpOut,im4);
        mobil2(atoms.coord,&lj,atoms.charge, inData.n_atoms,&mol, iic, &tmm[iic],&tmc[iic],&sdevpc,RanNum, fpOut,im2,&ifailc);
        #else
        mobil4(atoms.coord, &hs, &mol,RanNum,r,&ehsc[iic],&ehsm[iic],&pac[iic],&pam[iic],&imm,inData.n_atoms,fpOut,im4);
        mobil2(atoms.coord,&lj,atoms.charge, inData.n_atoms,&mol, iic, &tmm[iic],&tmc[iic],&sdevpc,RanNum, r,fpOut,im2,&ifailc);
        #endif

        if(imm<immmin) immmin=imm;
        if(imm>immmax) immmax=imm;  
    
        if (iic<inData.n_coord-1){
            ncoord(&inData,&mol,&atoms,&asympp[iic+1],fpIn,fpOut);
        } else fclose(fpIn); //nothing more to read

        //turn off the verbose output of mobil2/4
        im2=1;
        im4=1;
    }
    
    //if(inData.n_coord!=1) fclose(fpIn); //nothing more to read

    // Print Summary
    fprintf(fpOut,"\n\n\n SUMMARY\n\n");
    fprintf(fpOut," program version = 0.3 using %d threads\n",N_CORES);
    fprintf(fpOut," input file name = %s\n",inFile);
    fprintf(fpOut," input file label = %s\n",inData.label);
    if(tmmob!=0){
        if (inData.cUnit==EQUAL) fprintf(fpOut," using a uniform charge distribution\n\n");
        else if (inData.cUnit==CALC) fprintf(fpOut," using a calculated (non-uniform) charge distribution\n\n");
        else if (inData.cUnit==NONE) fprintf(fpOut," using no charge - only LJ interactions\n\n");
    }
    fprintf(fpOut," temperature =% -.4E\n",mol.temperature);

    #if TEST_CTYPES
    fprintf(fpOut," Getting Random Numbers from Python generator");
    #else
    fprintf(fpOut," using Mersenne Twister random number generator (GSL implementation:  %s ) with seed  = %ld\n",gsl_rng_name(r),Rseed);   
    #endif

    if (inData.n_coord==1){
        fprintf(fpOut," structural asymmetry parameter =% -.4E\n",asympp[0]);
        if(ehsm[0]!=0){ // NOTE : this line was in the original code, but there is no way that ehsm = 0 
            fprintf(fpOut,"\n mobility calculation by MOBIL4 (HS scattering)\n");
            fprintf(fpOut,"\n number of Monte Carlo trajectories =%7d\n",INUM*N_CORES);
            fprintf(fpOut," maximum number of reflections encountered =%3d\n",imm);
            fprintf(fpOut,"\n inverse average PA mobility =% -.4E\n",1.0/pam[0]);
            fprintf(fpOut," average PA cross section =% -.4E\n",pac[0]*1.e20);
            fprintf(fpOut,"\n inverse average EHS mobility =% -.4E\n",1.0/ehsm[0]);
            fprintf(fpOut," average EHS cross section =% -.4E\n",ehsc[0]*1e20);
        }
        
        if (tmm[0]!=0){ // NOTE : this line was in the original code, but there is no way that tmm = 0 
            fprintf(fpOut,"\n mobility calculation by MOBIL2 (trajectory method)\n");
            fprintf(fpOut,"\n trajectory parameters\n");
            fprintf(fpOut," sw1 =% -.4E       sw2 =% -.4E\n",SW1,SW2);
            fprintf(fpOut," dtsf1 =% -.4E     dtsf2 =% -.4E\n",DTSF1,DTSF2);
            fprintf(fpOut,"\n number of complete cycles (itn) =%6d\n",ITN*N_CORES);
            fprintf(fpOut," number of velocity points (inp) =%6d\n",INP);
            fprintf(fpOut," number of random points (imp) =%6d\n",IMP);
            fprintf(fpOut," total number of points =%7d\n",ITN*INP*IMP*N_CORES);
            fprintf(fpOut,"\n inverse average (second order) TM mobility =% -.4E\n",1.0/(tmm[0]));
            fprintf(fpOut," average TM cross section =% -.4E\n",tmc[0]*1.e20);
            fprintf(fpOut," standard deviation (percent) =% -.4E\n",sdevpc);
            fprintf(fpOut," number of failed trajectories =%4d\n",ifailc);
        }

} else{     // multiple coordinate sets
    if (ehsm[0]!=0){
        fprintf(fpOut,"\n mobility calculation by MOBIL4 (HS scattering)\n");
        fprintf(fpOut,"\n number of Monte Carlo trajectories =%7d\n",INUM*N_CORES);
        fprintf(fpOut,"\n minimum and maximum number of reflections =%3d  %3d\n",immmin,immmax);
        fprintf(fpOut,"\n\n    set     PA CS       PA MOB^-1      EHSS CS      EHSS MOB^-1    ASYMP\n");
        //DIFF: in the following printing the last element asympp is divided by 10 in the fortran code, but not here.
        // That is only a formatting hack, not needed here.
        for(i=0;i<inData.n_coord;i++) fprintf(fpOut," %5d   % -.4E   % -.4E   % -.4E   % -.4E    % -.4f\n",i+1,pac[i]*1.0e20,1.0/pam[i],ehsc[i]*1.0e20,1.0/ehsm[i],asympp[i]);
        
        for(i=0;i<inData.n_coord;i++) pacs+=pac[i]; 
        pacs/=(double) inData.n_coord;

        for(i=0;i<inData.n_coord;i++) pamob+=pam[i]; 
        pamob/=(double) inData.n_coord;
        
        for(i=0;i<inData.n_coord;i++) ehscs+=ehsc[i]; 
        ehscs/=(double) inData.n_coord;
        
        for(i=0;i<inData.n_coord;i++) ehsmob+=ehsm[i]; 
        ehsmob/=(double) inData.n_coord;
        
        for(i=0;i<inData.n_coord;i++) aasymp+=asympp[i]; 
        aasymp/=(double) inData.n_coord;
        
        //DIFF: in the following printing the last element aasymp is divided by 10 in the fortran code, but not here.
        // That is only a formatting hack, not needed here.
        fprintf(fpOut,"\n   AVGE  % -.4E   % -.4E   % -.4E   % -.4E    % -.4f\n",pacs*1.0e20,1.0/pamob,ehscs*1.0e20,1.0/ehsmob,aasymp);
    }
    if (tmm[0]!=0){
        fprintf(fpOut,"\n mobility calculation by MOBIL2 (trajectory method)\n");
        fprintf(fpOut,"\n trajectory parameters\n");
        fprintf(fpOut," sw1 =% -.4E       sw2 =% -.4E\n",SW1,SW2);
        fprintf(fpOut," dtsf1 =% -.4E     dtsf2 =% -.4E\n",DTSF1,DTSF2);
        fprintf(fpOut,"\n number of complete cycles (itn) =%6d\n",ITN*N_CORES);
        fprintf(fpOut," number of velocity points (inp) =%6d\n",INP);
        fprintf(fpOut," number of random points (imp) =%6d\n",IMP);
        fprintf(fpOut," total number of points =%7d\n",ITN*INP*IMP*N_CORES);
        fprintf(fpOut,"\n total number of failed trajectories =%4d\n",ifailc);
        fprintf(fpOut,"\n\n    set     TM CS       TM MOB^-1\n");
        for(i=0;i<inData.n_coord;i++) fprintf(fpOut," %5d   % -.4E   % -.4E\n",i+1,tmc[i]*1.0e20,1.0/tmm[i]);
        
        for(i=0;i<inData.n_coord;i++) tmcs+=tmc[i]; 
        tmcs/=(double) inData.n_coord;
        
        for(i=0;i<inData.n_coord;i++) tmmob+=tmm[i]; 
        tmmob/=(double) inData.n_coord;
    
        fprintf(fpOut,"\n   AVGE  % -.4E   % -.4E\n",tmcs*1.0e20,1.0/tmmob);
    }
}
    fclose(fpOut); //nothing more to write
    // Free some stuff
    FreeLj(&lj);
    FreeAtomData(&atoms);
    FreeHS(&hs);
    free(tmm);
    free(tmc);
    free(asympp);
    free(ehsc);
    free(ehsm);
    free(pac);
    free(pam);

    #if !TEST_CTYPES
    gsl_rng_free(r);
    #endif

    return(EXIT_SUCCESS);
}
