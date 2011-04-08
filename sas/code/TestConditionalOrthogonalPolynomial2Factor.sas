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
* Conditional power for multipvariate with two within subject factors
*/
TITLE1 "Test in multivariate model with two within factors - based on example 9 from POWERLIB";
TITLE2 "From Coffey C.S. and Muller K.E. (2003) Properties of ";
TITLE3 "internal pilots with the univariate approach to repeated measures";
TITLE4 "Statistics in Medicine, 22(15)";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&DATA_DIRECTORY";

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

ALPHA = .04;
OPT_ON = {UN HF GG BOX  HLT PBT WLK UMETHOD MMETHOD NOPRINT};
OPT_OFF = {COLLAPSE WARN};
ROUND = 15;

BETASCAL = 1;
THETA = {.25}#{.5 1 -1 .5}; * =Theta(cr) from 1st sentence *after* 
                            *  equation 7, Coffey and Muller (2003); 

* Following from Table II in Coffey and Muller (2003) *;
VARSTARE = {.47960 .01000 .01000 .01000}; * epsilon ~ .28 *; 
VARSTARF = {.34555 .06123 .05561 .04721}; * epsilon ~ .50 *;
VARSTARG = {.23555 .17123 .05561 .04721}; * epsilon ~ .72 *;
VARSTARH = {.12740 .12740 .12740 .12740}; * epsilon = 1 *;
VARSTAR = VARSTARE//VARSTARF//VARSTARG//VARSTARH;

SIGSCAL = {0.50 1.00 2.00}; * <=> gamma in Coffey and Muller (2003) *;

* Log base 2 spacing Clip (2,4,16) and Region(2,8,32) *;
* Get orthonormal U matrices *;
RUN UPOLY2({1 2 4},"A", {1 3 5},"B",
            UA,NMA, UB,NMB, UAB,NMAB);
U = UAB;
C = 1;

ESSENCEX = {1};
REPN = {20};

* Run with Muller Barton approximation ;
UCDF = J(4,1,1);
UMETHOD = J(2,1,1);
  DO IVAR = 1 TO 4 BY 1;
    SIGSTAR = DIAG(VARSTAR[IVAR,*]);
    SIGMA = U*SIGSTAR*U`;  * 1st paragraph in section 2.4, Coffey and Muller 2003 *;
    BETA = THETA*U`;       * 1st paragraph in section 2.4, Coffey and Muller 2003 *; 

    RUN POWER;
    HOLDALL = HOLDALL//_HOLDPOWER;
  END;
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF'};
filename out "&DATA_DIRECTORY\TestConditionalOrthogonalPolynomial2FactorMB.xml";
RUN powerResultsToXML(out, HOLDALL, TEST_LIST, 0);


* Run with Muller, Edwards, Simpson, Taylor approximation ;
FREE HOLDALL;
UCDF = J(4,1,2);
UMETHOD = J(2,1,2);
  DO IVAR = 1 TO 4 BY 1;
    SIGSTAR = DIAG(VARSTAR[IVAR,*]);
    SIGMA = U*SIGSTAR*U`;  * 1st paragraph in section 2.4, Coffey and Muller 2003 *;
    BETA = THETA*U`;       * 1st paragraph in section 2.4, Coffey and Muller 2003 *; 

    RUN POWER;
    HOLDALL = HOLDALL//_HOLDPOWER;
  END;
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF'};
filename out "&DATA_DIRECTORY\TestConditionalOrthogonalPolynomial2FactorMEST.xml";
RUN powerResultsToXML(out, HOLDALL, TEST_LIST, 0);


QUIT;


