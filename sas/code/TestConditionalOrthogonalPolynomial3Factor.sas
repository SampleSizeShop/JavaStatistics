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
* Conditional power for all tests with fixed design, triply repeated measures
*
* based on example 8 from POWERLIB (Johnson et al., 2009)
*/
TITLE1 "Conditional Power, Illustrate use of the UPOLY3 module";
TITLE2 "Factorial design, repeated measures: A, B, C between, D, E, F within";
%INCLUDE "common.sas";

LIBNAME DATA_DIR "&DATA_DIRECTORY";

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;
RESET FUZZ NOAUTONAME FW=6 LINESIZE=80;

ALPHA = .05;

* Choose dimensions of design *;
GA = 3; * =# groups for between factor A *;
GB = 3; * =# groups for between factor B *;
GC = 3; * =# groups for between factor C *;
TD = 3; * =#Times for within factor D *;
TE = 3; * =#Times for within factor E *;
TF = 3; * =#Times for within factor F *;

P = TD#TE#TF;
Q = GA#GB#GC;
ESSENCEX = I(Q);
BETA = J(Q,P,0);
BETA[1,1] = 1;
REPN = DO(2,12,2);
SIGMA = DIAG(DO(1,P,1)); * Variances are 1,2,3,...p *;

* Get orthonormal U matrices *;
CALL UPOLY3 ( (1:TD),"D", (1:TE),"E",  (1:TF),"F",
		          UD,UDLBL,   UE,UELBL,    UF,UFLBL, 
                 UDE,UDELBL, UDF,UDFLBL,  UEF,UEFLBL,  UDEF,UDEFLBL );
 
* Get orthonormal C matrices *;
CALL UPOLY3 ((1:GA),"A" , (1:GB),"B" , (1:GC),"C",
		         U1,CALBL,    U2,CBLBL,    U3,CCLBL,
                U12,CABLBL, U13,CACLBL,   U23,CBCLBL,  U123,CABCLBL);

BETASCAL = {9 18 27};
ROUND = 15;

OPT_ON = {UN HF GG BOX  HLT PBT WLK MMETHOD  UMETHOD MMETHOD NOPRINT};
OPT_OFF = {COLLAPSE};
BUG=" ";

C = U1`;
U = UD;
RUN POWER;
HOLDA=HOLDA//_HOLDPOWER;
PRINT / "AxD";
PRINT HOLDA[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

C = U12`;
U = UDE;
RUN POWER;
HOLDABDE = HOLDABDE//_HOLDPOWER;
PRINT / "AxB x DxE Interaction";
PRINT HOLDABDE[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

C = U123`;
U = UDEF;
RUN POWER;
HABCDEF = HABCDEF//_HOLDPOWER;
PRINT / "AxBxC x DxExF Interaction";
PRINT HABCDEF[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

HOLDALL = HOLDA//HOLDABDE//HABCDEF;
PRINT (NROW(HOLDALL));
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&DATA_DIRECTORY\TestConditionalOrthogonalPolynomial3Factor.xml";
RUN powerResultsToXML(out, HOLDALL, TEST_LIST, 0);

QUIT;


