/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2010 Regents of the University of Colorado.  
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package edu.cudenver.bios.power.test.general;

import java.text.DecimalFormat;
import java.util.List;

import edu.cudenver.bios.power.OneSampleStudentsTPower;
import edu.cudenver.bios.power.OneSampleStudentsTPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters;
import jsc.distributions.Normal;
import junit.framework.TestCase;

/**
 * Unit test for one sample students t-test
 * 
 */
public class TestPowerStudentT extends TestCase
{
	public static final String TITLE = "One Sample Student T";
    public static final int SIMULATION_SIZE = 10000;
    public static final boolean SIMULATE = true;
    public static final double PRECISION = 0.01;
    private static Normal normalDist = new Normal();
    private static DecimalFormat Number = new DecimalFormat("#0.0000");

	OneSampleStudentsTPowerParameters params = new OneSampleStudentsTPowerParameters();
	
	/**
	 * Build the parameters
	 */
    public void setUp()
    {
    	params.addAlpha(0.05);
    	params.addAlpha(0.01);
    	params.addMeans(0, 1);
    	params.addMeans(20, 24);
    	params.addVariance(4); 
    	params.addSampleSize(10);
    	params.addSampleSize(15);
    	params.addSampleSize(30);
    	params.addPower(0.1);
    	params.addPower(0.5);
    	params.addPower(0.8);
    	params.addPower(0.9);
    }
    
    /**
     * Test one tailed test against simulation
     */
    public void testPowerOneTailed()
    {
    	params.setTwoTailed(false);
    	
		System.out.println(TITLE + ": one tailed");
		int mismatches = checkPower(params);
		assertEquals(0, mismatches);
    }
    
    /**
     * Test two tailed test against simulation
     */
    public void testPowerTwoTailed()
    {
    	params.setTwoTailed(true);
		System.out.println(TITLE + ": two tailed");
		int mismatches = checkPower(params);
		assertEquals(0, mismatches);
    }
    
    /**
     * Check for errors during sample size run
     */
    public void testSampleSizeTwoTailed()
    {
    	params.setTwoTailed(true);
        OneSampleStudentsTPowerCalculator calc = new OneSampleStudentsTPowerCalculator();
        
        try
        {
        	System.out.println("Testing two tailed sample size");
        	System.out.println("Calculating sample sizes...");
        	List<Power> results = calc.getSampleSize(params);
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
    }
    
    
    /**
     * Run the power calculations for the specified parameters and tests
     * and assert whether they match simulation
     * 
     * @param params
     * @param testList
     * @returns the number of powers that failed to match simulation
     */
    public int checkPower(OneSampleStudentsTPowerParameters params)
    {
    	// create a power calculator
    	OneSampleStudentsTPowerCalculator calc = new OneSampleStudentsTPowerCalculator();
    	int simulationMatches = 0;

    	// perform the calculations
    	System.out.println("Calculating power...");
    	long startTime = System.currentTimeMillis();
    	List<Power> results = calc.getPower(params);
    	long calcTime = System.currentTimeMillis() - startTime;
    	System.out.println("Done.  Elapsed time: " +  ((double) calcTime / (double) 1000) + " seconds");
    	
    	System.out.println("Calc Power\tSim Power (p-value)\tTotal N\tDifference\tSigmat\tAlpha");

    	// perform the simulation if requested
    	List<Power> simResults = null;
    	if (SIMULATE)
    	{
    		System.out.println("Simulating power...");
        	startTime = System.currentTimeMillis();
    		simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
    		long simTime = System.currentTimeMillis() - startTime;
        	System.out.println("Done.  Elapsed time: " +  ((double) simTime / (double) 1000) + " seconds");
    	}
    	
    	// accumulate results
    	for(int i = 0; i < results.size(); i++)
    	{
    		OneSampleStudentsTPower power = (OneSampleStudentsTPower) results.get(i);
    		OneSampleStudentsTPower simPower = (simResults != null ? (OneSampleStudentsTPower) simResults.get(i) : null);
    		double simPvalue = Double.NaN;	

    		if (simPower != null)
    		{
    			simPvalue = zTest(power.getActualPower(), simPower.getActualPower());
    			if (simPvalue > PRECISION) simulationMatches++;
    		}
    		
    		System.out.println(Number.format(power.getActualPower()) + "\t" +
    				(simPower != null ? Number.format(simPower.getActualPower()) : "n/a")+ " (" + Number.format(simPvalue) + ")\t" +
    				Number.format(power.getTotalSampleSize()) + "\t" + 
    				Number.format(power.getDetectableDifference()) + "\t" + 
    				Number.format(power.getSigma()) + "\t" + 
    				Number.format(power.getAlpha()) + "\t");
    	}

    	int totalMismatches = 0;
    	if (SIMULATE)
    	{
    		totalMismatches += results.size() - simulationMatches;
    	}
    	return totalMismatches;
    }
    
    /**
     * compute p-value for comparing a calculated and simulated power by z-test 
     * 
     * @param calc calculated power value
     * @param sim simulated power value
     * @return p-value
     */
    private static double zTest(double calc, double sim)
    {
     	double z = Math.abs((sim - calc) / Math.sqrt((sim * (1 - sim)) / SIMULATION_SIZE));
        double p = 2 * normalDist.upperTailProb(z);
        return p;
    }
    
}
