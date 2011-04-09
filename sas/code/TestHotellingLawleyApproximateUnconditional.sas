/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2010 Regents of the University of Colorado.  
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/* 
* Test case for Hotelling Lawley, approximate, unconditional power
* Matches output in Table II, Glueck&Muller 2003
* Tests a multivariate GLH(F) for a design with
* 3 groups and a single baseline covariate
*/
TITLE "HLT approximate unconditional power";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&INPUT_DATA_DIRECTORY";

PROC IML SYMSIZE=2000 WORKSIZE=6000;
RESET SPACES=4;
/* Must include these files after call to PROC IML */
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEPA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEPE1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEQ1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\POWUA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\POWHA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\QLIB01.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEU4.IML"/NOSOURCE2;

USE DATA_DIR.BScales;
READ ALL 
     VAR _ALL_
     INTO   HOLDIN [ COLNAME =INNM ];
*INNM={"N" "PRQUANT" "TARGET" "FSCALE"};
CLOSE INOUT01.P0501;

DEBUG = 0;
IF DEBUG PRINT HOLDIN;

H_OR_U="H";
TEST="H";
P=4;
QF=3;
RHO=.50;
ALPHA=.05;

ZERO=1E-14;
ESSENCEF=I(QF);
CF=J(QF-1,1,-1)||I(QF-1);
A=NROW(CF);
CG=J(A,1,0);
C=CF||CG;
RANKX=NCOL(C);
U=I(P);
B=NCOL(U);
THETA0=J(A,B,0);
ESSBETAF=J(QF,P,0);
   ESSBETAF[1,1]=1;
   ESSBETAF[2,2]=2;
ESSBETAG=J(1,P-1,RHO)||{0};

IF DEBUG THEN PRINT C U, ESSBETAF ESSBETAG;

VARG=1; SIGMAG=VARG;
SIGMAY=I(P);
SIGMAE=SIGMAY-ESSBETAG`*SIGMAG*ESSBETAG;
SIGSTAR=U`*SIGMAE*U;
SIGSTVAL=EIGVAL(SIGSTAR);
IF DEBUG THEN PRINT SIGSTVAL;
IF MIN(SIGSTVAL)<ZERO THEN STOP;
ESSTD=C*(ESSBETAF//ESSBETAG)*U-THETA0;


SIGMAGY = SIGMAG#ESSBETAG;
SIGMAYG = SIGMAGY`;
IF DEBUG THEN PRINT SIGMAG, SIGMAY, SIGMAGY, SIGMAYG, ESSBETAG;

MAXPDIFF=.001;
MAXITER=50;
PROBLIST={0.5 };

DO ROW = 1 to NROW(HOLDIN);
  N      =HOLDIN[ROW,1];
  DELTA = HOLDIN[ROW,4];
  BETA = (DELTA#ESSBETAF)//ESSBETAG;
  IF DEBUG THEN PRINT BETA;
  THETA = C*BETA*U;
  THETAD=THETA-THETA0;
    IF DEBUG THEN PRINT THETA, THETA0;
  REPN=N/NCOL(ESSENCEF);
  DO IPROB=1 TO NCOL(PROBLIST);
	/* approximate unconditional power */
	APOWUNH1=BASEU4("A",H_OR_U,ESSENCEF,REPN,CF,CG,VARG,THETAD,SIGSTAR,
             ALPHA,MAXPDIFF);
    POWER=POWER//(DELTA||N||APOWUNH1[1,1]);
   END;
END;

NAMES = {'Beta-Scale' 'Total N' 'HLT Power'};
PRINT POWER[COLNAME=NAMES];

/* write to XML file */
filename out "&OUTPUT_DATA_DIRECTORY\TestHotellingLawleyApproximateUnconditional.xml"; 
file out;
	put "<powerList>";
	do i=1 to NROW(POWER);
		put "<glmmPower test='hlt' alpha='" @;
		put ALPHA @;
		put "' nominalPower='" @;
		put (POWER[i,3]) @;
		put "' actualPower='" @;
		put (POWER[i,3]) @;
		put "' betaScale='" @;
		put (POWER[i,1]) @;
		put "' sigmaScale='1' sampleSize='" @;
		put (POWER[i,2]) @;
		put "' powerMethod='unconditional' />";
	end;
	put "</powerList>";
closefile out;

QUIT;

