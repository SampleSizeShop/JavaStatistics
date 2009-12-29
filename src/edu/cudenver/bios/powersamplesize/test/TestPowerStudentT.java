/**
 * Unit test for one sample students t-test
 * 
 */
package edu.cudenver.bios.powersamplesize.test;

import java.util.ArrayList;

import edu.cudenver.bios.powersamplesize.PowerOneSampleStudentsT;
import edu.cudenver.bios.powersamplesize.parameters.SimplePowerSampleSizeParameters;
import junit.framework.TestCase;

public class TestPowerStudentT extends TestCase
{
    public static final int SIMULATION_SAMPLE_SIZE = 10000;
    public static final double PRECISION = 0.01;

    private ArrayList<SimplePowerSampleSizeParameters> testCaseParameters = null;
    
    /* list of mu0, muA, sigma, alpha, sampleSize for simulating power */
    private double[][] simulationTestCases = 
    { 
            {0,1,1,0.05,50}, 
            {0,-1,1,0.05,50}, 
            {20,22,4,0.05,34} 
    };

    public void setUp()
    {
        testCaseParameters = new ArrayList<SimplePowerSampleSizeParameters>();
        // build some test cases
        for(int i = 0; i < simulationTestCases.length; i++)
        {
            double[] simulation = simulationTestCases[i];
            SimplePowerSampleSizeParameters params = new SimplePowerSampleSizeParameters();
            params.setMu0(simulation[0]);
            params.setMuA(simulation[1]);
            params.setSigma(simulation[2]);
            params.setAlpha(simulation[3]);
            params.setSampleSize((int) simulation[4]);
            testCaseParameters.add(params);
        }
        
        
    }
    
    public void testPowerOneTailed()
    {
        int successCount = 0;
        PowerOneSampleStudentsT tpower = new PowerOneSampleStudentsT();
        
        try
        {
            for(SimplePowerSampleSizeParameters params: testCaseParameters)
            {
                params.setOneTailed(true);

                double simulatedPower = tpower.getSimulatedPower(params, SIMULATION_SAMPLE_SIZE);
                double power = tpower.getCalculatedPower(params);
                System.out.println("One-tail Simulated power: " + simulatedPower + ", Calculated power: " + power);
                if (Math.abs(simulatedPower - power) < PRECISION) successCount++;
            }

        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, simulationTestCases.length);
        }

    }
    
    public void testPowerTwoTailed()
    {
        int successCount = 0;
        PowerOneSampleStudentsT tpower = new PowerOneSampleStudentsT();
        
        try
        {
            for(SimplePowerSampleSizeParameters params: testCaseParameters)
            {
                params.setOneTailed(false);

                double simulatedPower = tpower.getSimulatedPower(params, SIMULATION_SAMPLE_SIZE);
                double power = tpower.getCalculatedPower(params);
                System.out.println("One-tail Simulated power: " + simulatedPower + ", Calculated power: " + power);
                if (Math.abs(simulatedPower - power) < PRECISION) successCount++;
            }

        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, simulationTestCases.length);
        }

    }
}
