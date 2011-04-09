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
* Conditional power for all tests with fixed design
*/
TITLE "Conditional Power, Fixed Design";
%INCLUDE "common.sas";

***************************************************************;
* Perform power calculations for a two sample T test,          ;
* replicating the results in "Increasing scientific power with ;
* statistical power", by K.E. Muller and V.A. Benignus,        ;
* Neurotoxicology and Teratology, vol 14, May-June, 1992       ;
* The code reports power for a limited number of predicted     ;
* differences in means, compared to the number of values       ;
* needed for plotting.                                         ;
* Based on example 1 from POWERLIB (Johnson et al, 2009)       ;
***************************************************************;

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

* Define inputs to power program;
ALPHA = 0.05;
SIGMA = {1};
SIGSCAL = {0.32 1.00 2.05};

ESSENCEX = I(2);
REPN = { 10 };

BETA = {0 1}`;
BETASCAL = DO(0,2.5,0.05);
C = {1 -1};

OPT_OFF= {C U};
ROUND = 15;
RUN POWER;

/* write the data to an XML file */
TEST_LIST = {"unirep"};
filename out "&OUTPUT_DATA_DIRECTORY\TestConditionalTwoSampleTTest.xml"; 
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 1);

QUIT;
