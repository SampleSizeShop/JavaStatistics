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
TITLE "Conditional Power, Fixed Design, Paired T-Test";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&DATA_DIRECTORY";
***************************************************************;
* Power for a paired t-test.  Based on example 2 from 
*   Johnson J.L., Muller K.E., Slaughter J.C., Gurka M.J., Gribbin M.J. and Simpson S.L. 
*   (2009) POWERLIB: SAS/IML software for computing power in multivariate linear models, 
*   Journal of Statistical Software, 30(5), 1-27.
***************************************************************;

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

* Define inputs to power program;
ALPHA = 0.05;
ESSENCEX = I(1);
SIGMA = {2 1, 1 2};
BETA = {0 1};
C = {1};
U = {1 -1}`;

SIGSCAL = {1};
BETASCAL = DO(0,2.5,0.5);
REPN = { 10 };
OPT_ON = {COLLAPSE};
OPT_OFF= {C U};

ROUND = 15;
RUN POWER;

/* write the data to an XML file */
TEST_LIST = {"unirep"};
filename out "&DATA_DIRECTORY\TestConditionalPairedTTest.xml"; 
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 1);

QUIT;
