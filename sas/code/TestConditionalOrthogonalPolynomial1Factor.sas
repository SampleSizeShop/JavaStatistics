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
* Orthogonal polynomial U matrix for 1 factor.  Based on example 7
* from POWERLIB (Johnson et al., 2009)
*/
TITLE "Orthogonal polynomial U matrix for 1 factor";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&DATA_DIRECTORY";

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

OPT_ON = {ORTHU UN HF GG BOX HLT PBT WLK MMETHOD UMETHOD MMETHOD};
OPT_OFF = {COLLAPSE};

ALPHA = .05;
VARIANCE = 1.5;  
RHO = 0.25;
* Create compound symmetric covariance structure *;
SIGMA = VARIANCE#(I(5)#(1-RHO) + J(5,5,RHO)); 

ESSENCEX = I(2);
REPN = {10,20,40};

BETASCAL = 1;
BETA = {0 0 0 0 1,
        1 0 0 0 0};
C = {1 -1};

TIMES ={2 4 6 8 10};
RUN UPOLY1(TIMES  ,"Time", USCORE , SCORENM );
U = USCORE;

* full precision;
ROUND = 15;
RUN POWER;

/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&DATA_DIRECTORY\TestConditionalOrthogonalPolynomial1Factor.xml";
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 0);

QUIT;
