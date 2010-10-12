TITLE1 "Confirming equality with JavaStats library";

***************************************************************;
* Perform power calculations for a two sample T test,          ;
* replicating the results in "Increasing scientific power with ;
* statistical power", by K.E. Muller and V.A. Benignus,        ;
* Neurotoxicology and Teratology, vol 14, May-June, 1992       ;
* The code reports power for a limited number of predicted     ;
* differences in means, compared to the number of values       ;
* needed for plotting.  Code for plot is in ExampleA02.sas     ;
***************************************************************;



PROC IML SYMSIZE=1000 WORKSIZE=2000;

%LET LMDIRECT = C:\Documents and Settings\kreidles\My Documents\Glimmpse\KeithMullerCode\LINMOD33\SOURCE\    ;
%INCLUDE "&LMDIRECT.LINMOD.IML" /   NOSOURCE2 ;
%INCLUDE "C:\Documents and Settings\kreidles\My Documents\Glimmpse\KeithMullerCode\powerlib\Iml\POWERLIB203.IML"/NOSOURCE2;
%INCLUDE "C:\Documents and Settings\kreidles\My Documents\Glimmpse\KeithMullerCode\kem\research\randomx\baseline\iml\SIMLIB1.IML";

SEEDE=3863;
SEEDG=1187;

/* This is the exact matrix produced on a single run with seed=3863
// to be cut/pasted into Java code for testing as needed
//        double[][] errorData = {
//                {-1.108957, 0.5362602, -0.875543},
//                {1.0065903, 1.2013761, 0.0718426},
//                { -0.88725, -0.197967,  0.293882},
//                {1.6828316, -0.164961, -0.275763},
//                {0.0914157, -1.490488, -1.231211},
//                {0.9088365, 0.8968534, 0.7452088},
//                {1.0151016, -0.892324, 1.2530911},
//                {-1.714659, -0.428179, -0.406727},
//                {-1.705252, 0.7142231, 0.0971129},
//                {0.4429702, 0.0308012, 0.1096798},
//                {-0.243068, 0.4415986,  0.362475},
//                {2.5378809, 0.8984633, 1.2051368},
//                { 1.3168095,  0.454879,  -0.49463},
//                { -1.274457, -0.644826, -1.013701},
//                { -1.13208, 1.7034088, 2.0165851},
//                {-0.001806, 0.7899111, -0.441356},
//                {  -0.192865, -1.040685, -0.476169},
//                {  -0.15294, 0.4456642, 0.5015633},
//                { -1.320593, -0.366892, 0.2203713},
//                { -0.349038, -0.791849, -0.425987},       
//        };
//        RealMatrix error = new Array2DRowRealMatrix(errorData);
//        
//        return error;
*/

P=3;
Q=4;
C=J(Q-1,1,1)||(-I(Q-1));
U=( J(P-1,1,1)||(-I(P-1)) )`;
*U={1 0 0, 0 1 0, 0 0 1};

ALPHA=.05;

VARIANCE=1;  RHO=0.4;
SIGMAE=VARIANCE#(I(P)#(1-RHO) + J(P,P,RHO)); *Compound symmetry;
*SIGMAE={1 0.2 0.3,0.2 1 0.2,0.3 0.2 1};
PRINT SIGMAE;
ESSENCEX=I(Q); *REPN={5, 10}; *REPN={5,10}; 
BETA=J(Q,P,0);
BETA[1,1]=2;
FSIGMAET=ROOT(SIGMAE);

OPT_ON = {"MULTTEST" "UNIRFORC" "NOPRINT"};
OPT_OFF={ "RSQUARED" "ECORR" "HEIVAL" "CANVEC" 
          "CANRSQ" "EVEC2" "EXTHETA" "MATTHETA" "UNITHETA" "MSH" "MSE"}; 
RUN SETOPT;
FREE OPT_OFF;

DEPVARS=NAMELIST("Y",1,P,1);
INDVARS=NAMELIST("f",1,Q,1);*||{"g"};
ZNAMES=DEPVARS||INDVARS;

N = 20;
MUG=J(N,1,0);
REPN={5};
F=ESSENCEX@J(REPN[1,1],1,1);
X=F;
XBETA=X*BETA;

* Simulation iterations ;
NREP = 10000;
* Rejection counts for unirep, gg, wilks, hotelling, pillai;
REJECT = {0 0 0 0 0};
DO REPLICAT=1 TO NREP;
    E=GAUSS1(N,J(N,P,0),FSIGMAET,SEEDE);
	PRINT E;
    Y=XBETA+E;
    Z=Y||X;
    RUN MAKESS;
    RUN FITMODEL;
    RUN TESTGLH;
	IF _URESUL_[1,5] <= alpha THEN REJECT[1,1]=REJECT[1,1]+1;
	IF _URESUL_[2,5] <= alpha THEN REJECT[1,2]=REJECT[1,2]+1;
	IF _STMAT1_[2,5] <= alpha THEN REJECT[1,3]=REJECT[1,3]+1;
	IF _STMAT1_[3,5] <= alpha THEN REJECT[1,4]=REJECT[1,4]+1;
	IF _STMAT1_[4,5] <= alpha THEN REJECT[1,5]=REJECT[1,5]+1;
END;
POWERHAT=(REJECT)/NREP;
PRINT  FSIGMAET;
PRINT _SIGMA_;
PRINT  _URESUL_;
PRINT _STMAT1_;
*HOLD=HOLD//(HOLDIN[ICASE,*]||DESIGN||POWERHAT);  

* calculate the powers for these tests with powerlib ;
OPT_OFF={ALPHA};
OPT_ON = {ORTHU UN GG WLK HLT PBT };
SIGMA = SIGMAE;
SIGSCAL={1};
BETA[1,1]=1;
BETASCAL=(2);
RUN POWER; 

PRINT "SIMULATED";
PRINT "UN	GG	WL	HLT	PBT";
PRINT POWERHAT;
QUIT;







