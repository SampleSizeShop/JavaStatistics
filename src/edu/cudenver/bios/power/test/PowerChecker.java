package edu.cudenver.bios.power.test;

import java.text.DecimalFormat;
import java.util.List;

import jsc.distributions.Normal;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class PowerChecker 
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double POWER_SIMULATION_COMPARISON_ALPHA = 0.01;
    private static Normal normalDist = new Normal();
    private static DecimalFormat Number = new DecimalFormat("#0.0000");
    
    /**
     * Run the power calculations for the specified parameters and tests
     * and assert whether they match simulation
     * 
     * @param params
     * @param testList
     * @returns the number of powers that failed to match simulation
     */
    public static int checkPower(GLMMPowerParameters params, 
            boolean univariate, boolean fixed)
    {
    	return PowerChecker.checkPower(params, univariate, fixed, true);   	
    }
    
    /**
     * Run the power calculations for the specified parameters and tests
     * and assert whether they match simulation
     * 
     * @param params
     * @param testList
     * @returns the number of powers that failed to match simulation
     */
    public static int checkPower(GLMMPowerParameters params, 
            boolean univariate, boolean fixed, boolean simulate)
    {
    	// create a power calculator
    	GLMMPowerCalculator calc = new GLMMPowerCalculator();

    	// build prefix string
    	String uniStr = (univariate ? "U" : "M");
    	String fixedStr = (fixed ? "F" : "R");

    	int matches = 0;

    	System.out.println("Calculating power...");
    	List<Power> results = calc.getPower(params);
    	System.out.println("Done.");
    	if (simulate)
    	{
    		System.out.println("Simulating power...");
    		List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
    		System.out.println("Done.");

    		System.out.println("P-value\tM/U\tF/R\tCalc Power\tSim Power\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");

    		for(int i = 0; i < results.size(); i++)
    		{
    			GLMMPower power = (GLMMPower) results.get(i);
    			GLMMPower simPower = (GLMMPower) simResults.get(i);
    			double pValue = zTest(power.getActualPower(), simPower.getActualPower());
    			if (pValue > POWER_SIMULATION_COMPARISON_ALPHA) matches++;
    			System.out.println(Number.format(pValue) + "\t" + uniStr + "\t" + fixedStr + "\t" + 
    					Number.format(power.getActualPower()) + "\t" +
    					Number.format(simPower.getActualPower()) + "\t" +  
    					power.getTest().toString() + "\t" + 
    					Number.format(power.getSigmaScale()) + "\t" + 
    					Number.format(power.getBetaScale()) + "\t" + 
    					power.getTotalSampleSize() + "\t" + 
    					Number.format(power.getAlpha()) + "\t" + 
    					power.getPowerMethod() + "\t" + 
    					power.getQuantile() + "\t");
    		}


    	}
    	else
    	{
    		System.out.println("M/U\tF/R\tCalc Power\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");
    		for(int i = 0; i < results.size(); i++)
    		{
    			GLMMPower power = (GLMMPower) results.get(i);
    			System.out.println(uniStr + "\t" + fixedStr + "\t" + 
    					Number.format(power.getActualPower()) + "\t" +  
    					power.getTest().toString() + "\t" + 
    					Number.format(power.getSigmaScale()) + "\t" + 
    					Number.format(power.getBetaScale()) + "\t" + 
    					power.getTotalSampleSize() + "\t" + 
    					Number.format(power.getAlpha()) + "\t" + 
    					power.getPowerMethod() + "\t" + 
    					power.getQuantile() + "\t");
    		}
    		matches = results.size();
    	}

    	return results.size() - matches;
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
