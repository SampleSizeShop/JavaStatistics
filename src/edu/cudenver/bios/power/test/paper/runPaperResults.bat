REM
REM This file runs all Junit tests for examples described in the GLIMMPSE paper.  
REM
java -cp "lib/JavaStatistics.jar;thirdparty/junit/4.7/lib/junit-4.7.jar;thirdparty/ApacheCommonsMath/2.1/lib/commons-math-2.1.jar;thirdparty/jsc/1.0/lib/jsc.jar" edu.cudenver.bios.power.test.paper.PaperResultsTestSuite > results.out

