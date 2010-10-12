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
package edu.cudenver.bios.power.test;

/**
 * Helper class which runs power calculations for the java unit tests
 * and compares against simulated values and SAS output
 * @author Sarah Kreidler
 */
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
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
    
    private ArrayList<Result> checkerResults = new ArrayList<Result>();
    private Timer timer = new Timer();
    
    private class Timer
    {
    	public long calculationMilliseconds;
    	public long simulationMilliseconds;

    	public Timer() 
    	{
    		reset();
    	};
    	
    	public void reset()
    	{
    		calculationMilliseconds = 0;
        	simulationMilliseconds = 0;
    	}
    	
    	public void addCalculationTime(long time)
    	{
    		calculationMilliseconds += time;
    	}
    	
    	public void addSimulationTime(long time)
    	{
    		simulationMilliseconds += time;
    	}
    }
    
    private class Result
    {
    	public double calculatedPower;
    	public double simulatedPower;
    	public double calcSimPvalue;
    	public double sasPower;
    	public double calcSasPvalue;
    	public GLMMPowerParameters.Test test;
    	public double sigmaScale;
    	public double betaScale;
    	public int totalSampleSize;
    	public double alpha;
    	public GLMMPowerParameters.PowerMethod powerMethod;
    	public double quantile;   	
    	
    	public Result(double calculatedPower,
    			double simulatedPower,
    			double calcSimPvalue,
    			double sasPower,
    			double calcSasPvalue,
    			GLMMPowerParameters.Test test,
    			double sigmaScale,
    			double betaScale,
    			int totalSampleSize,
    			double alpha,
    			GLMMPowerParameters.PowerMethod powerMethod,
    			double quantile)
    	{
    		this.calculatedPower = calculatedPower;
    		this.simulatedPower = simulatedPower;
    		this.calcSimPvalue = calcSimPvalue;
    		this.sasPower = sasPower;
    		this.calcSasPvalue = calcSasPvalue;
    		this.test = test;
    		this.sigmaScale = sigmaScale;
    		this.betaScale = betaScale;
    		this.totalSampleSize = totalSampleSize;
    		this.alpha = alpha;
    		this.powerMethod = powerMethod;
    		this.quantile = quantile;
    	}
    };
    
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
     * @return the number of powers that failed to match simulation
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
    	long calcTime = System.currentTimeMillis() - startTime;
    	System.out.println("Done.  Elapsed time: " +  ((double) calcTime / (double) 1000) + " seconds");
    	timer.addCalculationTime(calcTime);
    	
    	// perform the simulation if requested
    	List<Power> simResults = null;
    	if (simulate)
    	{
    		System.out.println("Simulating power...");
        	startTime = System.currentTimeMillis();
    		simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
    		long simTime = System.currentTimeMillis() - startTime;
        	System.out.println("Done.  Elapsed time: " +  ((double) simTime / (double) 1000) + " seconds");
        	timer.addSimulationTime(simTime);
    	}
    	
    	// accumulate results
    	for(int i = 0; i < results.size(); i++)
    	{
    		GLMMPower power = (GLMMPower) results.get(i);
    		GLMMPower simPower = (simResults != null ? (GLMMPower) simResults.get(i) : null);
    		GLMMPower sasPower = findMatchingSASPower(power);
    		double simPvalue = Double.NaN;	
    		double sasPvalue = Double.NaN;	

    		if (simPower != null)
    		{
    			simPvalue = zTest(power.getActualPower(), simPower.getActualPower());
    			if (simPvalue > POWER_SIMULATION_COMPARISON_ALPHA) simulationMatches++;
    		}
    		if (sasPower != null)
    		{
    			sasPvalue = zTest(power.getActualPower(), sasPower.getActualPower());
    			if (sasPvalue > POWER_SIMULATION_COMPARISON_ALPHA) sasMatches++;
    		}

    		checkerResults.add(new Result(power.getActualPower(),
    			(simPower != null ? simPower.getActualPower() : Double.NaN),
    			simPvalue,
    			(sasPower != null ? sasPower.getActualPower(): Double.NaN),
    			sasPvalue,
    			power.getTest(),
    			power.getSigmaScale(),
    			power.getBetaScale(),
    			power.getTotalSampleSize(),
    			power.getAlpha(),
    			power.getPowerMethod(),
    			power.getQuantile()));
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
    
    public void reset()
    {
    	checkerResults.clear();
    	timer.reset();
    }
    
    /**
     * Write the results to an HTML file
     * @param filename
     */
    public void outputResults(String title, String filename)
    {
    	// output the results
    	StringBuffer buffer = new StringBuffer();
    	
    	buffer.append("<html><head></head><body><h1>" + title + "</h1>");
    	
    	buffer.append("<table border='1' cellpadding='5'>");
    	buffer.append("<tr><td>Calculation time</td><td>" + ((double) timer.calculationMilliseconds / 1000)
    			+ "</td></tr>");
    	buffer.append("<tr><td>Simulation time</td><td>" + ((double) timer.simulationMilliseconds / 1000)
    			+ "</td></tr></table><p></p>");
    	
    	buffer.append("<table border='1' cellpadding='5'>");
    	buffer.append("<tr><th>Calc Power</th><th>Sim Power (p-value)</th><th>SAS Power (p-value)</th><th>Test</th><th>SigmaScale</th><th>BetaScale</th><th>Total N</th><th>Alpha</th><th>PowerMethod</th><th>Quantile</th></tr>");

    	for(Result result: checkerResults)
    	{
    		buffer.append("<tr><td>" + Number.format(result.calculatedPower) + "</td><td>" +  
    				Number.format(result.simulatedPower));
    		if (result.calcSimPvalue < POWER_SIMULATION_COMPARISON_ALPHA)
    			buffer.append(" <font color='red'>(" + Number.format(result.calcSimPvalue) + ")</font></td><td>");
    		else
    			buffer.append(" (" + Number.format(result.calcSimPvalue) + ")</td><td>");
    		buffer.append(Number.format(result.sasPower));
    		if (result.calcSasPvalue < POWER_SIMULATION_COMPARISON_ALPHA)
    			buffer.append(" <font color='red'>(" + Number.format(result.calcSasPvalue) + ")</font></td><td>");
    		else
    			buffer.append(" (" + Number.format(result.calcSasPvalue) + ")</td><td>");
    		buffer.append(result.test + "</td><td>" + 
    				Number.format(result.sigmaScale) + "</td><td>" + 
    				Number.format(result.betaScale) + "</td><td>" + 
    				result.totalSampleSize + "</td><td>" + 
    				Number.format(result.alpha) + "</td><td>" + 
    				result.powerMethod + "</td><td>" + 
    				result.quantile + "</td><td></tr>");
    	}
    	
    	buffer.append("</table></body></html>");

    	FileWriter writer = null;
    	try
    	{
    		writer = new FileWriter(filename);
    		writer.write(buffer.toString());
    	}
    	catch (Exception e)
    	{
    		// TODO:
    	}
    	finally
    	{
    		if (writer != null) try {writer.close(); } catch (Exception e) {};
    	}
    	
    }
    
    /**
     * Write the current result set to stdout.
     */
    public void outputResults()
    {
    	// output the results
    	System.out.println("Calculation time: " + ((double) timer.calculationMilliseconds / 1000));
    	System.out.println("Simulation time: " + ((double) timer.simulationMilliseconds / 1000));
    	System.out.println("Calc Power\tSim Power (p-value)\tSAS Power (p-value)\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");
    	
    	for(Result result: checkerResults)
    	{
    		System.out.println(Number.format(result.calculatedPower) + "\t" +  
    				Number.format(result.simulatedPower) + " (" + Number.format(result.calcSimPvalue) + ")\t" +
    				Number.format(result.sasPower) + " (" +  Number.format(result.calcSasPvalue) + ")\t" +
    				result.test + "\t" + 
    				Number.format(result.sigmaScale) + "\t" + 
    				Number.format(result.betaScale) + "\t" + 
    				result.totalSampleSize + "\t" + 
    				Number.format(result.alpha) + "\t" + 
    				result.powerMethod + "\t" + 
    				result.quantile + "\t");
    	}
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
