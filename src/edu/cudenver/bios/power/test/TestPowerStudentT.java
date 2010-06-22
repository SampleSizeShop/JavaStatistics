/**
 * Unit test for one sample students t-test
 * 
 */
package edu.cudenver.bios.power.test;

import java.util.List;

import edu.cudenver.bios.power.OneSampleStudentsTPower;
import edu.cudenver.bios.power.OneSampleStudentsTPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters;
import junit.framework.TestCase;

public class TestPowerStudentT extends TestCase
{
    public static final int SIMULATION_SAMPLE_SIZE = 10000;
    public static final double PRECISION = 0.01;
    
	OneSampleStudentsTPowerParameters params = new OneSampleStudentsTPowerParameters();

    public void setUp()
    {
    	params.addAlpha(0.05);
    	params.addAlpha(0.01);
    	params.addMeans(0, 1);
    	params.addMeans(20, 24);
    	params.addVariance(5); 
    	params.addSampleSize(10);
    	params.addSampleSize(15);
    	params.addSampleSize(100);
    	params.addPower(0.1);
    	params.addPower(0.5);
    	params.addPower(0.8);
    	params.addPower(0.9);
    }
    
    public void testPowerOneTailed()
    {
    	params.setOneTailed(true);
        int successCount = 0;
        OneSampleStudentsTPowerCalculator calc = new OneSampleStudentsTPowerCalculator();
        
        int numResults = 0;
        try
        {
        	System.out.println("Testing one tailed power");
        	System.out.println("Calculating power values...");
        	List<Power> results = calc.getPower(params);
        	System.out.println("Simulation power values...");
        	List<Power> simulationResults = calc.getSimulatedPower(params, SIMULATION_SAMPLE_SIZE);
        	numResults = results.size();
        	for(int i = 0; i < results.size() && i < simulationResults.size(); i++)
        	{
        		OneSampleStudentsTPower calcResult = (OneSampleStudentsTPower) results.get(i);
        		OneSampleStudentsTPower simulationResult = (OneSampleStudentsTPower) simulationResults.get(i);
        		System.out.println("CALCULATED: " + calcResult.toXML());
        		System.out.println("SIMULATED: " + simulationResult.toXML());
                if (Math.abs(simulationResult.getPower() - calcResult.getPower()) < PRECISION) successCount++;
        	}
        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, numResults);
        }
    }
    
    public void testPowerTwoTailed()
    {
    	params.setOneTailed(false);
        int successCount = 0;
        OneSampleStudentsTPowerCalculator calc = new OneSampleStudentsTPowerCalculator();
        
        int numResults = 0;
        try
        {
        	System.out.println("Testing two tailed power");
        	System.out.println("Calculating power values...");
        	List<Power> results = calc.getPower(params);
        	System.out.println("Simulation power values...");
        	List<Power> simulationResults = calc.getSimulatedPower(params, SIMULATION_SAMPLE_SIZE);
        	numResults = results.size();
        	for(int i = 0; i < results.size() && i < simulationResults.size(); i++)
        	{
        		OneSampleStudentsTPower calcResult = (OneSampleStudentsTPower) results.get(i);
        		OneSampleStudentsTPower simulationResult = (OneSampleStudentsTPower) simulationResults.get(i);
        		System.out.println("CALCULATED: " + calcResult.toXML());
        		System.out.println("SIMULATED: " + simulationResult.toXML());
                if (Math.abs(simulationResult.getPower() - calcResult.getPower()) < PRECISION) successCount++;
        	}
        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            assertEquals(successCount, numResults);
        }
    }
    
    public void testSampleSizeTwoTailed()
    {
    	params.setOneTailed(false);
        int successCount = 0;
        OneSampleStudentsTPowerCalculator calc = new OneSampleStudentsTPowerCalculator();
        
        int numResults = 0;
        try
        {
        	System.out.println("Testing two tailed sample size");
        	System.out.println("Calculating sample sizes...");
        	List<Power> results = calc.getSampleSize(params);
        	numResults = results.size();
        	for(int i = 0; i < results.size(); i++)
        	{
        		OneSampleStudentsTPower calcResult = (OneSampleStudentsTPower) results.get(i);
        		System.out.println("CALCULATED: " + calcResult.toXML());
        	}
        }
        catch(Exception e)
        {
            System.out.println(e.getMessage());
            fail();
        }
        finally
        {
            //assertEquals(successCount, numResults);
        }
    }
}
