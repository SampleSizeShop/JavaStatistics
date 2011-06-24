@echo off
REM
REM This file runs all Junit tests for examples described in the GLIMMPSE paper.  
REM
echo =====================;
echo Starting replication of Table 2. 
echo Please note that simulations from this table may take several hours to complete.
java -cp "lib/JavaStatistics.jar;thirdparty/junit/4.7/lib/junit-4.7.jar;thirdparty/ApacheCommonsMath/2.1/lib/commons-math-2.1.jar;thirdparty/jsc/1.0/lib/jsc.jar" edu.cudenver.bios.power.test.paper.PaperTable2TestSuite > tables2.out
echo Done replicating results from Table 2.  Please see tables2.out for results.
echo =====================;
echo Starting replication of Table 3. 
echo Please note that simulations from this table may take several hours to complete.
java -cp "lib/JavaStatistics.jar;thirdparty/junit/4.7/lib/junit-4.7.jar;thirdparty/ApacheCommonsMath/2.1/lib/commons-math-2.1.jar;thirdparty/jsc/1.0/lib/jsc.jar" edu.cudenver.bios.power.test.paper.PaperTable3TestSuite > tables3.out
echo Done replicating results from Table 3.  Please see tables3.out for results.
echo =====================;
echo Starting replication of power results from Figure 4.
java -cp "lib/JavaStatistics.jar;thirdparty/junit/4.7/lib/junit-4.7.jar;thirdparty/ApacheCommonsMath/2.1/lib/commons-math-2.1.jar;thirdparty/jsc/1.0/lib/jsc.jar" org.junit.runner.JUnitCore edu.cudenver.bios.power.test.paper.TestLawANOVA > figure4.out
echo Done replicating results from Figure 4.  See figure4.out for results.
echo =====================;
echo Starting replication of power results from Figure 5.
java -cp "lib/JavaStatistics.jar;thirdparty/junit/4.7/lib/junit-4.7.jar;thirdparty/ApacheCommonsMath/2.1/lib/commons-math-2.1.jar;thirdparty/jsc/1.0/lib/jsc.jar" org.junit.runner.JUnitCore edu.cudenver.bios.power.test.paper.TestLawANCOVA > figure5.out
echo Done replicating results from Figure 5.  See figure5.out for results.>>>>>>> .merge-right.r929
