/* Test case for non-centrality distribution - for comparison
* with the Java Statistics library
*/
TITLE "Test: Davies Algorithm for CDF of Weighted Sum of non-central Chi-Squares";
%INCLUDE "common.sas";

PROC IML SYMSIZE=2000 WORKSIZE=6000;
RESET SPACES=4;
%INCLUDE "io.iml" / NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEPA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEPE1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEQ1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\POWUA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\POWHA1.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\QLIB01.IML"/NOSOURCE2;
%INCLUDE "&GLUECK_MULLER_IML_DIRECTORY\BASEU4.IML"/NOSOURCE2;

LAMBDAS = {7,-3,5};
OMEGAS = {10,2,1};
DFS = {1,2,1};
QUANTILES = {10 20 30 40 50 60 70 80 90 100};
ACCURACY = 0.001;

do I = 1 to NCOL(QUANTILES);
	Q = QUANTILES[I];
	P = QPROB(Q, LAMBDAS, DFS, OMEGAS, ACCURACY);
	HOLDPROB = HOLDPROB//(Q||P);
end;

NAMES = {"Quantile" "Pr(Q < q)"};
PRINT HOLDPROB[COLNAME=NAMES];


/* Write the input chi square information and probability results to an xml file */
filename out "&DATA_DIRECTORY\TestDaviesAlgorithm.xml"; 
file out;
	put "<testCase accuracy='" @; 
	put ACCURACY @;
	put "'>";
	put "<chiSquareList>";
	do i=1 to NROW(LAMBDAS);
		put "<chiSquare lambda='" @;
		put (LAMBDAS[i]) @;
		put "' omega='" @;
		put (OMEGAS[i]) @;
		put "' df='" @;
		put (DFS[i]) @;
		put "' />";
	end;
	put "</chiSquareList>";
	put "<cdfResultList>";
	do i=1 to NROW(HOLDPROB);
		put "<cdfResult quantile='" @;
		put (HOLDPROB[i,1]) @;
		put "' probability='" @;
		put (HOLDPROB[i,2]) @;
		put "' />";
	end;
	put "</cdfResultList>";
	put "</testCase>";
closefile out;

QUIT;


