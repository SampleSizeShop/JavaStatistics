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


FREE HOLDALL;
DO CLTYPE = 1 to 2;
  RUN CIPOWER;
END;
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestUnivariateConfidenceIntervals.xml"; 
RUN powerAndCIResultsToXML(out, HOLDALL, TEST_LIST, 1);


/*
* Repeat for a multivariate design: paired t-test, beta known
*/
* Define inputs to power program (only resetting values that changed from univariate);
ESSENCEX = I(1);
SIGMA = {2 1, 1 2};
BETA = {0 1};
C = {1};
U = {1 -1}`;

OPT_ON={UN BOX HF GG WLK PBT HLT ALPHA CLTYPE BETASCAL SIGSCAL ALPHA_CL ALPHA_CU NOPRINT};
OPT_OFF={COLLAPSE WARN};
ROUND = 15;

FREE HOLDALL;
CLTYPE = 1;
RUN CIPOWER;
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestMultivariateConfidenceIntervalsBetaKnown.xml"; 
RUN powerAndCIResultsToXML(out, HOLDALL, TEST_LIST, 0);


/*
* Repeat for a multivariate design: paired t-test, beta unknown
*/
FREE HOLDALL;

* Define inputs to power program (only resetting values that changed from univariate;
ESSENCEX = I(1);
SIGMA = {2 1, 1 2};
BETA = {0 1};
C = {1};
U = {1 -1}`;

OPT_ON={WLK PBT HLT ALPHA CLTYPE BETASCAL SIGSCAL ALPHA_CL ALPHA_CU NOPRINT};
OPT_OFF={COLLAPSE WARN};

ROUND = 15;
RUN POWER;
FREE HOLDALL;
CLTYPE = 2;
RUN CIPOWER;
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write to XML file */
TEST_LIST = {'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestMultivariateConfidenceIntervalsBetaUnknown.xml"; 
RUN powerAndCIResultsToXML(out, HOLDALL, TEST_LIST, 0);

QUIT;
