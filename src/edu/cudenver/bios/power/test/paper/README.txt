#
# Instructions for reproducing results and examples from the
# GLIMMPSE manuscript.  This script replicates the following results:
#
# - GLMM(F) accuracy results from Table 2
# - GLMM(F,g) accuracy results from Table 3
#
# Note that simulation values and CPU times may vary depending 
# on the processor and operating system.
#

-------------------------------------------------------
System requirements

- Java Runtime Environment 1.7.0_x
- Windows 7

Note, the examples may run on other OS platforms, but have only been tested in Windows 7.

--------------------------------------------------------
Instructions for running examples

1. Download the GLIMMPSEresults.zip from 

http://samplesizeshop.org/files/2012/12/GLIMMPSEresults.zip

2. Unzip the package.  By default, this should create a GLIMMPSEresults directory
in the same location as the GLIMMPSEresults.zip file.

3. Change the working directory to the GLIMMPSEresults directory

4. Run the runPaperResults.bat file

On Windows 7, this will open a DOS window to run the script.  It will 
produce 4 output files:
- table2.out  : statistics displayed in Table 2
- table3.out  : statistics displayed in Table 3

Please note that the simulations used for Tables 2 and 3 are time consuming,
so the script will likely take anywhere from 12-24 hours to complete depending
on processor speed.

------------------------------------------------------
Third party libraries

All required thirdparty libraries are included in the GLIMMPSEresults.zip
distribution.  License information for the supporting libraries is included
with each library in the thirdparty directory.

-----------------------------------------------------

Please contact Sarah Kreidler if you encounter any problems running
this script.

Author: Sarah Kreidler
Email: sarah.kreidler@ucdenver.edu
