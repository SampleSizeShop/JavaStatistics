package edu.cudenver.bios.power.test;

import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.xml.sax.InputSource;

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

    private boolean simulate = true;
    private List<GLMMPower> sasPowers = null;
    
    public PowerChecker() {}
    
    public PowerChecker(boolean compareAgainstSimulation)
    {
    	this.simulate = compareAgainstSimulation;
    }
    
    public PowerChecker(String sasOutputFile, boolean compareAgainstSimulation)
    throws IllegalArgumentException
    {
    	this.simulate = compareAgainstSimulation;
    	FileReader reader = null;
    	// parse the sas xml file
    	try
    	{
            DocumentBuilderFactory factory =  DocumentBuilderFactory.newInstance();

            DocumentBuilder builder = factory.newDocumentBuilder();
            reader = new FileReader(sasOutputFile);
            Document doc = builder.parse(new InputSource(reader));
            sasPowers = SASOutputParser.parsePowerResults(doc);
    	}
    	catch (Exception e)
    	{
    		throw new IllegalArgumentException("Parsing of SAS XML failed: " + e.getMessage());
    	}
    	finally
    	{
    		if (reader != null) try {reader.close();} catch (Exception e) {};
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
    public int checkPower(GLMMPowerParameters params)
    {
    	// create a power calculator
    	GLMMPowerCalculator calc = new GLMMPowerCalculator();
    	int simulationMatches = 0;
    	int sasMatches = 0;

    	// perform the calculations
    	System.out.println("Calculating power...");
    	long startTime = System.currentTimeMillis();
    	List<Power> results = calc.getPower(params);
    	long endTime = System.currentTimeMillis();
    	System.out.println("Done.  Elapsed time: " +  ((double) (endTime - startTime) / (double) 1000) + " seconds");
    	
    	// perform the simulation if requested
    	List<Power> simResults = null;
    	if (simulate)
    	{
    		System.out.println("Simulating power...");
        	startTime = System.currentTimeMillis();
    		simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
        	endTime = System.currentTimeMillis();
        	System.out.println("Done.  Elapsed time: " +  ((double) (endTime - startTime) / (double) 1000) + " seconds");
    	}
    	
    	// output the results
    	System.out.println("Calc Power\tSim Power\tP-value\tSAS Power\tP-value\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");
    	
    	for(int i = 0; i < results.size(); i++)
    	{
    		GLMMPower power = (GLMMPower) results.get(i);
    		GLMMPower simPower = (simResults != null ? (GLMMPower) simResults.get(i) : null);
    		GLMMPower sasPower = findMatchingSASPower(power);
    			
    		System.out.print(Number.format(power.getActualPower()) + "\t");
    		if (simPower != null)
    		{
    			double pvalue = zTest(power.getActualPower(), simPower.getActualPower());
    			if (pvalue > POWER_SIMULATION_COMPARISON_ALPHA) simulationMatches++;
    			System.out.print(Number.format(simPower.getActualPower()) + "\t" +  Number.format(pvalue) + "\t");
    		}
    		else
    		{
    			System.out.print("n/a\tn/a\t");
    		}
    		
    		if (sasPower != null)
    		{
    			double pvalue = zTest(power.getActualPower(), sasPower.getActualPower());
    			if (pvalue > POWER_SIMULATION_COMPARISON_ALPHA) sasMatches++;
    			System.out.print(Number.format(sasPower.getActualPower()) + "\t" +  Number.format(pvalue) + "\t");
    		}
    		else
    		{
    			System.out.print("n/a\tn/a\t");
    		}
    		
    		System.out.println(power.getTest().toString() + "\t" + 
    				Number.format(power.getSigmaScale()) + "\t" + 
    				Number.format(power.getBetaScale()) + "\t" + 
    				power.getTotalSampleSize() + "\t" + 
    				Number.format(power.getAlpha()) + "\t" + 
    				power.getPowerMethod() + "\t" + 
    				power.getQuantile() + "\t");
    	}

    	int totalMismatches = 0;
    	if (simulate)
    	{
    		totalMismatches += results.size() - simulationMatches;
    	}
    	if (sasPowers != null)
    	{
    		totalMismatches += results.size() - sasMatches;
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
    
    private GLMMPower findMatchingSASPower(GLMMPower power)
    {
    	if (sasPowers != null)
    	{
    		GLMMPower match = null;
    		for(GLMMPower sasPower: sasPowers)
    		{
    			if (power.getTest() == sasPower.getTest()  &&
    					power.getAlpha() == sasPower.getAlpha() &&
    					power.getTotalSampleSize() == sasPower.getTotalSampleSize() &&
    					power.getBetaScale() == sasPower.getBetaScale() &&
    					power.getSigmaScale() == sasPower.getSigmaScale() &&
    					power.getPowerMethod() == sasPower.getPowerMethod() &&
    					(Double.isNaN(sasPower.getQuantile()) || power.getQuantile() == sasPower.getQuantile()))
    			{
    				match = sasPower;
    				break;
    			}
    		}
    		return match;
    	}
    	else
    	{
    		return null;
    	}
    }
}
