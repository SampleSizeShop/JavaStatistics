/**
 * Unit test for one sample students t-test
 * 
 */
package edu.cudenver.bios.powersamplesize.test;

import edu.cudenver.bios.powersamplesize.PowerOneSampleStudentsT;
import edu.cudenver.bios.powersamplesize.SampleSizeOneSampleStudentsT;
import edu.cudenver.bios.powersamplesize.parameters.SimplePowerSampleSizeParameters;
import junit.framework.TestCase;

public class TestSampleSizeStudentT extends TestCase
{

    /* list of mu0, muA, sigma, alpha, sampleSize for simulating power */
    private double[][] sampleSizeTestCases = 
    { 
            {0,1,1,0.05,0.8}, 
            {0,-1,1,0.05,0.6}, 
            {20,22,4,0.05,0.8} 
    };

    public void testSampleSizeOneTailed()
    {
        int successCount = 0;
        SampleSizeOneSampleStudentsT calc = new SampleSizeOneSampleStudentsT();
        
        try
        {
            for(int i = 0; i < sampleSizeTestCases.length; i++)
            {
                double[] simulation = sampleSizeTestCases[i];
                SimplePowerSampleSizeParameters params = new SimplePowerSampleSizeParameters();
                params.setMu0(simulation[0]);
                params.setMuA(simulation[1]);
                params.setSigma(simulation[2]);
                params.setAlpha(simulation[3]);
                params.setPower(simulation[4]);
                params.setOneTailed(true);
                
                int n = calc.getSampleSize(params);
                PowerOneSampleStudentsT powerCalc = new PowerOneSampleStudentsT();
                params.setSampleSize(n);
                System.out.println("One-tail sample size: " + n + " actual power: " + powerCalc.getCalculatedPower(params));
                if (n > 0) successCount++;
                
            }

        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, sampleSizeTestCases.length);
        }

    }
    
    public void testSampleSizeTwoTailed()
    {
        int successCount = 0;
        SampleSizeOneSampleStudentsT calc = new SampleSizeOneSampleStudentsT();
        
        try
        {
            for(int i = 0; i < sampleSizeTestCases.length; i++)
            {
                double[] simulation = sampleSizeTestCases[i];
                SimplePowerSampleSizeParameters params = new SimplePowerSampleSizeParameters();
                params.setMu0(simulation[0]);
                params.setMuA(simulation[1]);
                params.setSigma(simulation[2]);
                params.setAlpha(simulation[3]);
                params.setPower(simulation[4]);
                params.setOneTailed(false);
                
                int n = calc.getSampleSize(params);
                PowerOneSampleStudentsT powerCalc = new PowerOneSampleStudentsT();
                params.setSampleSize(n);
                System.out.println("Two-tail sample size: " + n + " actual power: " + powerCalc.getCalculatedPower(params));
                if (n > 0) successCount++;
            }

        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, sampleSizeTestCases.length);
        }
    }
}
