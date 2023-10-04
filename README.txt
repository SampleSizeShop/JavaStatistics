  Java Statistics.  A java library providing power/sample size estimation for 
 the general linear model.
  
 Copyright (C) 2010 Regents of the University of Colorado.  
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
------------------------------
1. INTRODUCTION
------------------------------

This library provides power and sample size calculations for the
general linear multivariate model.  It is a component of the
Glimmpse software system (http://glimmpse.samplesizeshop.org/)

The power calculations are based on the work of Professor Keith E. Muller
and colleagues.  A full list of related publications are available at:

http://samplesizeshop.org/education/related-publications/

------------------------------
2.  LATEST VERSION
------------------------------

Version 2.0.0

------------------------------
3.  DOCUMENTATION
------------------------------

Documentation is available from the project web site:

http://samplesizeshop.org/documentation/glimmpse/

------------------------------
4. DEPENDENCIES
------------------------------

Java Runtime Environment 1.7.0 or higher
JUnit 4.7
Apache Commons Math 3.0 or higher
JSC Statistics Package (http://www.jsc.nildram.co.uk/)
Apache Ant 1.8.1

------------------------------
5.  SUPPORT
------------------------------

This library is provided without warranty.

For questions regarding this library, please email sarah.kreidler@ucdenver.edu

------------------------------
6.  ANT BUILD SCRIPT
------------------------------

The main build.xml script is located in the ${JAVASTATISTICS_HOME}/build
directory.  To compile the library, cd to ${JAVASTATISTICS_HOME}/build
and type

cd to ${JAVASTATISTICS_HOME}/build
ant

The resulting jar file is called

${JAVASTATISTICS_HOME}/lib/edu.ucdenver.bios.javastatistics-${version}.jar

The build script assumes that the a directory called thirdparty is
installed at the same directory level as ${JAVASTATISTICS_HOME}.

------------------------------
7. TEST / DEMO PROGRAMS
------------------------------

Several unit test / example programs are provided under the
${JAVASTATISTICS_HOME}/src directory in the package

edu.cudenver.bios.power.test.published

There are associated SAS (requires SAS 9.2 or higher) programs in the 
${JAVASTATISTICS_HOME}/sas/code directory which run the equivalent
SAS/IML implementation of the test cases.  The SAS test cases produce
XML in ${JAVASTATISTICS_HOME}/sas/data which brew install are then read by the
Java unit tests.

The following test cases are available:

* TestConditionalMultivariate.java - conditional power for multivariate inputs with fixed effects only
* TestConditionalUnivariate.java - conditional power for univariate inputs with fixed effects only
* TestHotellingLawleyApproximateQuantile.java - approximate quantile power for the 
   Hotelling Lawley Trace for a multivariate design with a baseline covariate
* TestHotellingLawleyExactQuantile.java - exact (Davies' algorithm) quantile power for the 
   Hotelling Lawley Trace for a multivariate design with a baseline covariate
* TestHotellingLawleyApproximateUnconditional.java - approximate unconditional power for the 
   Hotelling Lawley Trace for a multivariate design with a baseline covariate
* TestHotellingLawleyExactUnconditional.java - exact (Davies' algorithm) unconditional power 
   for the Hotelling Lawley Trace for a multivariate design with a baseline covariate
* TestUnirepApproximateQuantile.java - approximate quantile power for the 
   Univariate Approach to Repeated Measures for a multivariate design with a baseline covariate
* TestUnirepExactQuantile.java - exact (Davies' algorithm) quantile power for the 
   Univariate Approach to Repeated Measures for a multivariate design with a baseline covariate
* TestUnirepApproximateUnconditional.java - approximate unconditional power for the 
   Univariate Approach to Repeated Measures for a multivariate design with a baseline covariate
* TestUnirepExactUnconditional.java - exact (Davies' algorithm) unconditional power 
   for the Univariate Approach to Repeated Measures for a multivariate design with a baseline covariate

Note that all tests of the "univariate approach to repeated measures" include uncorrected
results as well as Geisser-Greenhouse, Hyunh-Feldt, and Box corrected results.

------------------------------
8. KNOWN ISSUES
------------------------------

- The UNIREP tests for quantile and unconditional use the Muller, Edwards, Taylor 2004
approximation for the expected value of the unirep epsilon.  Therefore, the results are
closer to simulated values, but may not match the original SAS/IML implementation
based on Muller, Barton 1989.

- The calculated values for the one sample t-test match other software 
implementations (R, Russ Length's Power Calculator), but do not currently
match simulation.  There may be a problem in the t-test simulator.

------------------------------
9. CONTRIBUTORS / ACKNOWLEDGEMENTS
------------------------------

This library was created by Dr. Sarah Kreidler and Dr. Deb Glueck
at the University of Colorado Denver.

Special thanks to the following individuals were instrumental in completion of this project:
Dr. Keith Muller
Dr. Anis Karimpour-Fard
Dr. Jackie Johnson

