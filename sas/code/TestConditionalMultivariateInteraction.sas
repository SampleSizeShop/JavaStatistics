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
* Multivariate design with interaction test.  Based on 
* Example 5 from POWERLIB (Johnson et al. 2009)
*/
TITLE "Multivariate design with interaction hypothesis";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&INPUT_DATA_DIRECTORY";

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

OPT_ON = {ORTHU UN HF GG BOX HLT PBT WLK MMETHOD UMETHOD MMETHOD};
* Specifying the option ORTHU in OPT_ON allows the program to provide;
* an orthonormal U matrix if one is not given by the user;
* This is the case for the following code;

* Create contrast matrices C and U (non-orthonormal) *;
P = 3;
Q = 4;
C = J(Q-1,1,1)||(-I(Q-1));
U = ( J(P-1,1,1)||(-I(P-1)) )`;

ALPHA = .01;

VARIANCE = 1;  
RHO = 0.4;
SIGMA = VARIANCE#(I(P)#(1-RHO) + J(P,P,RHO)); *Compound symmetry;
SIGSCAL = {1, 2};

ESSENCEX = I(Q); 
REPN = {5,10};
BETA = J(Q,P,0);
BETA[1,1] = 1;
BETASCAL = DO(0, 2.0 , 0.50);

/* MMETHOD = {4,4,4}; * Two moment null approximations + OBrien and Shieh (1992) 
					 noncentrality multiplier ON;

UCDF = {4,2,2,4};  * UN and Box (4):
					 Exact via Davies' algorithm (1980), as in Muller, 
					 Edwards, Simpson, and Taylor (2007). If exact fails, 
					 then switch to approximation 2, MEST (2007);
                   * HF and GG (2): 
					 Muller, Edwards, Simpson, and Taylor (2007) 
					 approximation;
*/

* Output full precision ;
ROUND = 15;

RUN POWER;

/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestConditionalMultivariateInteraction.xml";
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 0);

QUIT;
