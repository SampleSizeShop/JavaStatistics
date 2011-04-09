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
* Conditional power
*/
TITLE "Conditional Power, Two Sample T-test with 3D plot";
%INCLUDE "common.sas";

***************************************************************;
* Perform power calculations for a two sample T test,          ;
* with a 3D plot showing trade-offs between                    ;
* based on EXAMPLE3.sas from powerlib                          ;
***************************************************************;

/* delete the power data set if it exists */
PROC DATASETS LIBRARY=WORK;
DELETE PWRDT1;
RUN; QUIT;

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

* Define inputs to power program;
OPT_ON = {DS};
* Output not printed to screen since NOPRINT option ON *;

ROUND = 15;
ALPHA = {.01};
SIGMA = {.068};  
SIGSCAL = {1} ;
BETA = {0 1}`;
C = {-1 1};
ESSENCEX = I(2);
REPN = DO(3,18,3);
BETASCAL = DO(0,.75,.05);
RUN POWER;

/* write the data to an XML file */
TEST_LIST = {"unirep"};
filename out "&OUTPUT_DATA_DIRECTORY\TestConditionalTwoSampleTTest3DPlot.xml"; 
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 1);

QUIT;

*** Section that creates plot - this section taken directly from EXAMPLE3.sas 
from powerlib 2.1 by Johnson et al. ***;

FILENAME OUT01 "&RESULTS_DIRECTORY\TestConditionalTwoSampleTTest3DPlot.png";

%MACRO GRAPH(TILT,ROTATE);
TITLE1;
PROC G3D DATA=ONE GOUT=OUT01;
PLOT BETASCAL*TOTAL_N=POWER/
     ZMIN=0 ZMAX=1.0 ZTICKNUM=6   YTICKNUM=4    XTICKNUM=6   SIDE;
LABEL TOTAL_N ="N"    BETASCAL="Delta"      POWER   ="Power" ;
RUN;
%MEND;

PROC G3GRID DATA=PWRDT1 OUT=ONE;
      GRID  BETASCAL*TOTAL_N=POWER /SPLINE NAXIS1=16 NAXIS2=11 ;

GOPTIONS GSFNAME=OUT01 DEVICE=PNG 
CBACK=WHITE COLORS=(BLACK) HORIGIN=0IN VORIGIN=0IN
HSIZE=5IN VSIZE=3IN HTEXT=12PT FTEXT=TRIPLEX;

%GRAPH(90, 300);

