package edu.cudenver.bios.powersamplesize.test;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.powersamplesize.PowerGLMM;
import edu.cudenver.bios.powersamplesize.SampleSizeGLMM;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;
import junit.framework.TestCase;

public class TestSampleSizeGLMM extends TestCase
{
    private static final double PRECISION = 0.01;
    private static final double ALPHA = 0.05;
    private static final double POWER = 0.80;
    
    private static final String[] testNames = {"None", "Unirep", "W","PB", "HLT"}; 
    
    public void testValidUnivariateMatrix()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs(true);

        PowerGLMM powerCalc = new PowerGLMM();
        SampleSizeGLMM calc = new SampleSizeGLMM();
        int tests = 0; // number of tests run
        int goodResults = 0; // number of matches between calculated and simulated power
        for(TestStatistic statistic: TestStatistic.values())
        {
            if (statistic == TestStatistic.NONE) continue; // skip the "NONE" statistic for now
            tests++;
            System.out.println("Good inputs, Test: " + testNames[tests]);

            goodParams.setTestStatistic(statistic);
            try
            {
                int sampleSize = (int) Math.ceil(calc.getSampleSize(goodParams));
                // get the actual power for that sample size
                RealMatrix design = goodParams.getDesignEssence().getFullDesignMatrix(sampleSize);
                goodParams.setDesign(design);
                goodParams.setSampleSize(sampleSize);
                double calcPower = powerCalc.getCalculatedPower(goodParams);
                // report results
                System.out.println("Sample size: " + sampleSize + " actual power=" + calcPower);
                if (Math.abs(calcPower - goodParams.getPower()) > PRECISION)
                    System.out.println("Unable to calculate power within tolerance of " + PRECISION);
                if (sampleSize > 0) goodResults++;
            }
            catch(Exception e)
            {
                System.out.println("Failed to calculate power: " + e.getMessage());
            }        
        }
        assertEquals(tests, goodResults);

    }
    
    
    private LinearModelPowerSampleSizeParameters buildValidUnivariateInputs(boolean essenceOnly)
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        params.setPower(POWER);
        
        // build beta matrix
        double [][] beta = {{0},{1.5}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{2.05}};
        params.setSigma(new Array2DRowRealMatrix(sigma));
        // build design matrix
        double[][] essenceData = {{1,0},{0,1}};
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(1,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(1,3));
        params.setDesignEssence(essenceMatrix);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
        
        return params;     
        
    }
}
