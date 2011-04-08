REM
REM This file runs all SAS code which generates XML inputs to the java unit tests.  Note, this only includes
REM examples described in the GLIMMPSE paper.  Additional .sas files must be run individuals.
REM
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalMultivariate.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalMultivariateInteraction.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalMultivariateWithConfidenceLimits.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalOrthogonalPolynomial1Factor.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalOrthogonalPolynomial2Factor.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalOrthogonalPolynomial3Factor.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalPairedTTest.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalTwoSampleTTest.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalTwoSampleTTest3DPlot.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalUnivariate.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestConditionalUnivariateWithConfidenceLimits.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestDaviesAlgorithm.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestHotellingLawleyApproximateQuantile.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestHotellingLawleyApproximateUnconditional.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestHotellingLawleyExactQuantile.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestHotellingLawleyExactUnconditional.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestUnirepApproximateQuantile.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestUnirepApproximateUnconditional.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestUnirepExactQuantile.sas
"C:\Program Files\SAS\SASFoundation\9.2\sas.exe" -NOSPLASH -NOLOG -NOPRINT -sysin ./TestUnirepExactUnconditional.sas
