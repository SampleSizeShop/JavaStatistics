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
TITLE "Test for univariate and multivariate confidence intervals";
%INCLUDE "common.sas";

***************************************************************;
* Confidence limits for a univariate model test
* based on example 4 from POWERLIB v21
* Creates Figure 1 in Taylor and Muller, 1995, Amer Statistician, 49, p43-47
***************************************************************;

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

* convenience routine to run power for various CI widths;
START CIPOWER;
  * Statements to create two-sided confidence limits *;
  N_EST = 24;        *# Obs for variance estimate;
  RANK_EST = 2;      *# model df for study giving variance estimate;
  ALPHA_CL = .025;   *Lower confidence limit tail size;
  ALPHA_CU = .025;   *Upper confidence limit tail size;
  ROUND = 15;
  RUN POWER;
  HOLDALL = HOLDALL // _HOLDPOWER;

  * Statements to create one-sided lower confidence limits *;
  ALPHA_CL = 0.05;   *Lower confidence limit tail size;
  ALPHA_CU =   0;    *Upper confidence limit tail size;
  *Since WORK.PWRDT1 already exists, WORK.PWRDT2 is created.;

  RUN POWER;
  HOLDALL = HOLDALL // _HOLDPOWER;

  * Statements to create one-sided upper confidence limits *;

  ALPHA_CL = 0;      *Lower confidence limit tail size;
  ALPHA_CU = 0.05;   *Upper confidence limit tail size;

  RUN POWER;
  HOLDALL = HOLDALL // _HOLDPOWER;
FINISH CIPOWER;

/*
* Univariate for both beta known and unknown
*/
* Define inputs to power program;
ESSENCEX = I(2);  *Balanced two group t test, cell mean coding;
REPN = {12};

BETA = {0 1}`;
BETASCAL = {0 0.25 0.5 0.75}; 

ALPHA={.01}; 
C    ={1 -1}; 

SIGMA  = {.068}; 
SIGSCAL = {1}   ; 
OPT_ON = {NOPRINT};
ROUND = 15;

FREE HOLDALL;
DO CLTYPE = 1 to 2;
  RUN CIPOWER;
END;
PRINT 'Univariate, includes beta known and unknown';
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestUnivariateConfidenceIntervals.xml"; 
RUN powerAndCIResultsToXML(out, HOLDALL, TEST_LIST, 1);


/*
* Repeat for a multivariate design, beta known
*/
* Define inputs to power program (only resetting values that changed from univariate);
ESSENCEX = I(3);
RHO = 0.4;
SIGMA = I(3)*(1-rho) + J(3,3,rho);
BETA = {1 0 0,
		0 0 0,
		0 0 0};
BETASCAL = {0.75 1 2};
C = {1 -1 0, 
	 1 0 -1};
U = C`;

OPT_ON={UN BOX HF GG WLK PBT HLT ALPHA CLTYPE BETASCAL SIGSCAL ALPHA_CL ALPHA_CU NOPRINT ORTHU};
OPT_OFF={COLLAPSE WARN};

FREE HOLDALL;
CLTYPE = 1;
RUN CIPOWER;
PRINT 'Multivariate, Beta Known';
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestMultivariateConfidenceIntervalsBetaKnown.xml"; 
RUN powerAndCIResultsToXML(out, HOLDALL, TEST_LIST, 0);


/*
* Repeat for a multivariate design, beta unknown
*/
FREE HOLDALL;
* disable UNIREP tests for CI's with beta and sigma estimated;
OPT_ON={WLK PBT HLT ALPHA CLTYPE BETASCAL SIGSCAL ALPHA_CL ALPHA_CU NOPRINT ORTHU};
OPT_OFF={UN BOX HF GG COLLAPSE WARN};

CLTYPE = 2;
RUN CIPOWER;
PRINT 'Multivariate, Beta Unknown';
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'wl' 'pbt' 'hlt'};
TEST_POWERCOL = {15,12,9};
filename out "&OUTPUT_DATA_DIRECTORY\TestMultivariateConfidenceIntervalsBetaUnknown.xml"; 
* had to hard code the XML here since you can't calculate CI's for;
* UNIREP tests when beta is unknown - the multivariate XML functions;
* used in the previous tests cases expect columns for all UNIREP/MULTIREP tests;
file out;
put "<powerList>";
do t = 1 to NCOL(TEST_LIST);
	POWERCOL = TEST_POWERCOL[t];
	do i=1 to NROW(HOLDALL);
		CLUPPERCOL = POWERCOL + 1;
		CLLOWERCOL = POWERCOL - 1;
		put "<glmmPower test='" @;
		put (TEST_LIST[t]) @;
		put "' alpha='" @;
		put (HOLDALL[i,1]) @;
		put "' nominalPower='" @;
		put (HOLDALL[i,POWERCOL]) best17. @;
		put "' actualPower='" @;
		put (HOLDALL[i,POWERCOL]) best17. @;
		put "' alphaLower='" @;
		put (HOLDALL[i,6]) @;
		put "' alphaUpper='" @;
		put (HOLDALL[i,7]) @;
		put "' ciLower='" @;
		put (HOLDALL[i,CLLOWERCOL]) best17. @;
		put "' ciUpper='" @;
		put (HOLDALL[i,CLUPPERCOL]) best17. @;
		put "' betaScale='" @;
		put (HOLDALL[i,3]) @;
		put "' sigmaScale='" @;
		put (HOLDALL[i,2]) @;
		put "' sampleSize='" @;
		put (HOLDALL[i,4]) @;
		put "' powerMethod='conditional' />";
    end;
end;
put "</powerList>";
closefile out;


QUIT;
