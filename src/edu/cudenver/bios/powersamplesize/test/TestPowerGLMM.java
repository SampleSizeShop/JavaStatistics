package edu.cudenver.bios.powersamplesize.test;

import org.apache.commons.math.linear.Array2DRowRealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.PowerGLMM;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;

import junit.framework.TestCase;

public class TestPowerGLMM extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double PRECISION = 0.01;
    private static final double ALPHA = 0.05;    
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 4.3;
    private static final double[] BETA_SCALE = {0, 2.5, 0.5};
    private static final double[] SIGMA_SCALE = {2};
    private static final int[] SAMPLE_SIZE = {10, 20};
    
    public void testValidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();

        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPower("Valid Univariate, Fixed, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Univariate, Fixed, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Univariate, Fixed, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Univariate, Fixed, W",calc, goodParams);        
    }
    
    public void testInvalidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();
        goodParams.setBeta(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Univariate, Fixed, UNIREP", calc, goodParams);
  
    }
    
    public void testValidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs();

        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPower("Valid Univariate, Fixed, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Univariate, Fixed, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Univariate, Fixed, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Univariate, Fixed, W",calc, goodParams);        
        
    }
    
    public void testInvalidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs();

        goodParams.setSigmaError(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Univariate, Fixed, UNIREP", calc, goodParams);
        
    }
    
    public void testValidUnivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();     
        
    }
    
    public void testInvalidUnivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();


        
    }
    
    public void testValidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs();

        ColumnMetaData randColMD = new ColumnMetaData();
        randColMD.setPredictorType(PredictorType.RANDOM);
        randColMD.setMean(MEAN);
        randColMD.setVariance(VARIANCE);
        goodParams.getDesignEssence().setColumnMetaData(2, new ColumnMetaData());
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPower("Valid Univariate, Fixed, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Univariate, Fixed, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Univariate, Fixed, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Univariate, Fixed, W",calc, goodParams);  


        
    }
    
    public void testInvalidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();


        
    }
    
    private void checkPowerFail(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        try
        {
            double calculated = calc.getCalculatedPower(params);
            double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) params, SIMULATION_SIZE);
            assert(Math.abs(simulated - calculated) < PRECISION);
        }
        catch(Exception e)
        {
            assert(true);
        }
            
    }
    
    private void checkPower(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        int tests = 0; // number of tests run
        int matches = 0; // number of matches between calculated and simulated power

        for(int sampleSize: SAMPLE_SIZE)
        {
            for(double sigmaScale : SIGMA_SCALE)
            {
                for(double betaScale : BETA_SCALE)
                {
                    tests++;     
                    LinearModelPowerSampleSizeParameters testParams = 
                        new LinearModelPowerSampleSizeParameters(params);
                    // scale the inputs
                    testParams.setBeta(params.getBeta().scalarMultiply(betaScale));
                    testParams.setSigmaError(params.getSigmaError().scalarMultiply(sigmaScale));
                    testParams.setSampleSize(sampleSize);
                    try
                    {
                        double calculated = calc.getCalculatedPower(testParams);
                        double simulated = calc.getSimulatedPower(testParams, SIMULATION_SIZE);

                        System.out.println(label + "["+betaScale+","+sigmaScale+","+sampleSize
                                +"]>> Calculated power: " + calculated + ", simulated power: " + simulated);
                        if (Math.abs(simulated - calculated) < PRECISION) matches++;
                    }
                    catch(Exception e)
                    {
                        System.out.println(label + "Failed to calculate power: " + e.getMessage());
                    }        

                }
            }
        }
        assertEquals(tests, matches);
    }
    
    private LinearModelPowerSampleSizeParameters buildValidUnivariateInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // build design matrix
        double[][] essenceData = {{1,0},{0,1}};
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(10,1));
        params.setDesignEssence(essenceMatrix);

        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
        
        return params;     
    }   

    
    private LinearModelPowerSampleSizeParameters buildValidMultivariateInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        double [][] beta = {{0},{1}, {1,2}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // build design matrix
        double[][] essenceData = {{1,0,0},{0,1,0},{0,0,1}};
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(2, new RowMetaData(10,1));
        params.setDesignEssence(essenceMatrix);

        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
        
        return params;     
    }   
}
