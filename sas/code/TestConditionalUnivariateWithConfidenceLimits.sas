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
TITLE "Conditional Power, 2 group t-test with confidence limits";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&DATA_DIRECTORY";
***************************************************************;
* Confidence limits for a univariate model test
* based on example 4 from POWERLIB v21
* Creates Figure 1 in Taylor and Muller, 1995, Amer Statistician, 49, p43-47
***************************************************************;

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

* Define inputs to power program;
ESSENCEX = I(2);  *Balanced two group t test, cell mean coding;
REPN = {12};

BETA = {0 1}`;
BETASCAL = DO(0,.75,.01); 

ALPHA={.01}; 
C    ={1 -1}; 

SIGMA  = {.068}; 
SIGSCAL = {1}   ; 

* Statements to create two-sided confidence limits *;

CLTYPE = 1;
N_EST = 24;        *# Obs for variance estimate;
RANK_EST = 2;      *# model df for study giving variance estimate;
ALPHA_CL = .025;   *Lower confidence limit tail size;
ALPHA_CU = .025;   *Upper confidence limit tail size;
ROUND = 15;
RUN POWER;

/* write to XML file */
TEST_LIST = {"unirep"};
filename out "&DATA_DIRECTORY\TestConditionalUnivariateWithConfidenceLimits.xml"; 
RUN powerAndCIResultsToXML(out, _HOLDPOWER, TEST_LIST, 1);


QUIT;
