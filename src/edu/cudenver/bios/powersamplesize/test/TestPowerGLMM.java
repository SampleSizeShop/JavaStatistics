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
    private static final int SIMULATION_SAMPLE_SIZE = 10000;
    private static final double PRECISION = 0.01;
    private static final double ALPHA = 0.05;
    private static final int SAMPLE_SIZE = 50;
    
    private static final String[] testNames = {"None", "Unirep", "W","PB", "HLT"}; 
    
    public void testValidUnivariateMatrix()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs(false);

        PowerGLMM calc = new PowerGLMM();
        int tests = 0; // number of tests run
        int matches = 0; // number of matches between calculated and simulated power
        for(TestStatistic statistic: TestStatistic.values())
        {
            if (statistic == TestStatistic.NONE) continue; // skip the "NONE" statistic for now
            tests++;
            System.out.println("Good inputs, Test: " + testNames[tests]);
            
            goodParams.setTestStatistic(statistic);
            try
            {
                double calculated = calc.getCalculatedPower(goodParams);
                double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) goodParams, SIMULATION_SAMPLE_SIZE);

                System.out.println("Calculated power: " + calculated + ", simulated power: " + simulated);
                if (Math.abs(simulated - calculated) < PRECISION) matches++;
            }
            catch(Exception e)
            {
                System.out.println("Failed to calculate power: " + e.getMessage());
            }        
        }
        assertEquals(tests, matches);
        
    }
    
    public void testValidUnivariateEssenceMatrix()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs(true);

        PowerGLMM calc = new PowerGLMM();
        int tests = 0; // number of tests run
        int matches = 0; // number of matches between calculated and simulated power
        for(TestStatistic statistic: TestStatistic.values())
        {
            if (statistic == TestStatistic.NONE) continue; // skip the "NONE" statistic for now
            tests++;
            System.out.println("Good inputs, Test: " + testNames[tests]);
            
            goodParams.setTestStatistic(statistic);
            try
            {
                double calculated = calc.getCalculatedPower(goodParams);
                double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) goodParams, SIMULATION_SAMPLE_SIZE);

                System.out.println("Calculated power: " + calculated + ", simulated power: " + simulated);
                if (Math.abs(simulated - calculated) < PRECISION) matches++;
            }
            catch(Exception e)
            {
                System.out.println("Failed to calculate power: " + e.getMessage());
            }        
        }
        assertEquals(tests, matches);
        
    }
    
    public void testNonPositiveDefiniteDesign()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildNonPositiveDefiniteInputs();
        PowerGLMM calc = new PowerGLMM();
        try
        {
            double calculated = calc.getCalculatedPower(goodParams);
            double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) goodParams, SIMULATION_SAMPLE_SIZE);
            
            assertEquals(calculated, simulated, PRECISION);
        }
        catch(Exception e)
        {
            fail();
        }   
    }
    
    public void testNonConformingMatrix()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildNonConformingInputs();
        PowerGLMM calc = new PowerGLMM();
        try
        {
            double calculated = calc.getCalculatedPower(goodParams);
            double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) goodParams, SIMULATION_SAMPLE_SIZE);
            
            assertEquals(calculated, simulated, PRECISION);
        }
        catch(Exception e)
        {
            fail();
        }   
    }
    
    private LinearModelPowerSampleSizeParameters buildValidUnivariateInputs(boolean essenceOnly)
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        double [][] beta = {{0},{2.0}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{2.05}};
        params.setSigma(new Array2DRowRealMatrix(sigma));
        // build design matrix
        if (!essenceOnly)
        {
            double [][] design = {
                    {1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},
                    {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}};
            params.setDesign(new Array2DRowRealMatrix(design));
        }
        else
        {
            double[][] essenceData = {{1,0},{0,1}};
            EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
            essenceMatrix.setRowMetaData(0, new RowMetaData(6,1));
            essenceMatrix.setRowMetaData(1, new RowMetaData(18,1));
//            ColumnMetaData cmd = new ColumnMetaData();
//            cmd.setMean(5);
//            cmd.setVariance(3);
//            cmd.setPredictorType(PredictorType.RANDOM);
//            essenceMatrix.setColumnMetaData(1, cmd);
            params.setDesignEssence(essenceMatrix);
        }
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
        
        return params;     
        
    }
    
    private LinearModelPowerSampleSizeParameters buildNonConformingInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        
        // build theta matrix
        
        // build sigma matrix
        
        //
        
        return params;
    }
    
    
    private LinearModelPowerSampleSizeParameters buildNonPositiveDefiniteInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        
        // build theta matrix
        
        // build sigma matrix
        
        //
        
        return params;
        
    }
    
    private LinearModelPowerSampleSizeParameters buildMissingInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        params.setSampleSize(SAMPLE_SIZE);
        
        return params;
    }
}
