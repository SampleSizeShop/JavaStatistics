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
* Conditional power for multivariate design with confidence limits 
* based on example 6 from POWERLIB (Johnson et al., 2009)
*
* Modifications from the original example:
* - we produce confidence limits for total n = 10,20,30,...,100
* - we produce results for all available tests
* - only a single data set is produced which is used to produce all of the plots
* - the sigma.sd2 data set is copied (and renamed) from p0104.sd2 file distributed
*   with POWERLIB
*/
TITLE "Conditional Power, Fixed Design, Multivariate with confidence limits";
%INCLUDE "common.sas";

LIBNAME IN V612 "&DATA_DIRECTORY\";

PROC DATASETS LIBRARY=WORK;
DELETE powerCL;
RUN; 

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

/* set output options - note, unlike example 6, we calculate power
* for all of the tests so the holdpower matrix will have appropriate
* dimensions for use with the XMLUtilities.iml module
*/
OPT_ON={UN BOX HF GG WLK PBT HLT ALPHA CLTYPE BETASCAL SIGSCAL ALPHA_CL ALPHA_CU NOPRINT};
OPT_OFF={COLLAPSE WARN};

* create power inputs;

* read in the covariance matrix from  ;
USE IN.sigma;
READ ALL VAR {ANT LEFT POST RIGHT} WHERE(_TYPE_ = "COV") INTO INSIGMA;

* Rounding and viewing covariance matrix *;
RNM={ANT LEFT POST RIGHT};
ALPHA = .05/6;
SIGMA = ROUND(INSIGMA,.0001);
PRINT SIGMA[FORMAT=11.5 COLNAME=RNM];

ESSENCEX = I(10);
REPN =DO(2,10,1);

P = 4;
Q = 2;
* Pattern of means for Gender by Region *;
BETARG= J(2,4, 3.2) + {.30}#({ -1  0 1  0,
                               -1  0 1  0}) ;	
PRINT BETARG[COLNAME=RNM];
BETASCAL={1};

* contrasts;
C = {1 -1} @ J(1,5,1);

REGION = {1,2,3,4};
RUN UPOLY1(REGION,"REGION",U1,REGU);
U=U1;

PRINT C, U;

* Output full precision ;
ROUND = 15;

* Describe dataset from whence the estimates came;
* needed to construct CI's *;

CLTYPE = 1;     *Estimated variance only, with fixed means*;
N_EST = 21;     *# Obs for variance estimate*;
RANK_EST = 1;   *# Model DF for study giving variance estimate*;
ALPHA_CL = .025;   	*Lower confidence limit tail size*;
ALPHA_CU = .025;   	*Upper confidence limit tail size*;

FREE HOLDALL;
DO DELTA=0 TO .20 BY .0008;
  * Creation of Beta matrix based on varying Gender differences *;
  BETARGD=BETARG + (J(2,2,0)||(DELTA//(-DELTA))||J(2,1,0));	
  * Final Beta matrix with age groups added *;
  BETA= BETARGD @ J(5,1,1) ;
  RUN POWER;
  HOLDALL=HOLDALL//( _HOLDPOWER||J(NROW(_HOLDPOWER),1,DELTA) );
END;

NAMES = _HOLDPOWERLBL || "DELTA";
CREATE powerCL FROM HOLDALL [COLNAME = NAMES]; 
APPEND FROM HOLDALL;
CLOSE powerCL;


PRINT HOLDALL[COLNAME=NAMES];
/* write the data to an XML file 
* note - we compute all of the powers to keep holdpower aligned
* but only output GG 
*/
TEST_LIST = {'unirepGG'};
filename out "&DATA_DIRECTORY\TestConditionalMultivariateWithConfidenceLimits.xml";
RUN powerAndCIResultsToXML(out, holdall, TEST_LIST, 0);
QUIT;

*** Section that creates plot with confidence limits for power for N=40***;

FILENAME OUT01 "&DATA_DIRECTORY\TestConditionalMultivariateWithConfidenceLimits.png";

GOPTIONS RESET=ALL GSFNAME=OUT01 DEVICE=PNG
CBACK=WHITE COLORS=(BLACK) HORIGIN=0IN VORIGIN=0IN
HSIZE=5IN VSIZE=3IN HTEXT=12PT FTEXT=TRIPLEX;

SYMBOL1 L=1 I=SPLINE V=NONE W=2;
SYMBOL2 L=20 I=SPLINE V=NONE W=2 R=2;


AXIS1 ORDER=(0 TO 1.0 BY .20)    W=3 MINOR=NONE MAJOR=(W=2)
      LABEL=(ANGLE=-90 ROTATE=90);
AXIS2 LABEL=("Delta") ORDER=(0 TO .20 BY .05) W=3 MINOR=NONE MAJOR=(W=2);

PROC GPLOT DATA=powerCL;
WHERE TOTAL_N=40;
PLOT (POWER_GG POWER_GG_L POWER_GG_U )*DELTA/ OVERLAY NOFRAME
       VZERO VAXIS=AXIS1 HZERO HAXIS=AXIS2 NOLEGEND;
LABEL POWER_GG ="Power" 
	  POWER_GG_L ="Power" 
	  POWER_GG_U ="Power";
RUN; 

