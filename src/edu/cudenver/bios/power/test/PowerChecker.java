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

import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.w3c.dom.Document;
import org.xml.sax.InputSource;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.utils.ConfidenceInterval;

public class PowerChecker 
{
	private class SASPower extends GLMMPower
	{
		double delta = Double.NaN;
		public SASPower(Test test, double alpha, 
				double nominalPower, double actualPower, int sampleSize,
				double betaScale, double sigmaScale, 
				GLMMPowerParameters.PowerMethod method,
				double quantile, ConfidenceInterval confidenceInterval, double delta)
		{
			super(test, alpha, nominalPower, actualPower, sampleSize,
					betaScale, sigmaScale, method, quantile, confidenceInterval);
			this.delta = delta;
		}
	}
	
	private double positivityThreshold = CholeskyDecompositionImpl.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
	private double symmetryThreshold = CholeskyDecompositionImpl.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;
    private static final int SIMULATION_SIZE = 10000;
    private static DecimalFormat Number = new DecimalFormat("#0.0000000");
    private static DecimalFormat LongNumber = new DecimalFormat("#0.00000000");

    private boolean simulate = true;
    private List<GLMMPower> sasPowers = null;
    
    // maximum allowed deviation
	private double sasTolerance = 0.00001;
	private double simTolerance = 0.05;
	
    // results
    private ArrayList<Result> checkerResults = new ArrayList<Result>();
    private double maxSasDeviation = 0;
    private double maxSimDeviation = 0;
    private double maxSaslowerCIDeviation = 0;
    private double maxSasUpperCIDeviation = 0;
    private int sasResultsIndex = 0;
    
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
    	GLMMPower calculatedPower;
    	double sasPower;
    	double sasCILower;
    	double sasCIUpper;
    	double simulatedPower;
    	double sasDeviation;
    	double sasCILowerDeviation;
    	double sasCIUpperDeviation;
    	double simulationDeviation;
    	
    	public Result(GLMMPower calculatedPower,
    			double sasPower,
    			double sasCILower,
    			double sasCIUpper,
    			double sasDeviation,
    			double sasCILowerDeviation,
    			double sasCIUpperDeviation,
    			double simulatedPower,
    			double simulationDeviation)
    	{
    		this.calculatedPower = calculatedPower;
    		this.sasPower = sasPower;
        	this.sasCILower = sasCILower;
        	this.sasCIUpper = sasCIUpper;
        	this.sasDeviation = sasDeviation;
        	this.sasCILowerDeviation = sasCILowerDeviation;
        	this.sasCIUpperDeviation = sasCIUpperDeviation;
    		this.simulatedPower = simulatedPower;
        	this.simulationDeviation = simulationDeviation;
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
    public void checkPower(GLMMPowerParameters params)
    {
    	// create a power calculator
    	GLMMPowerCalculator calc = new GLMMPowerCalculator();
    	calc.setPositivityThreshold(positivityThreshold);
    	calc.setSymmetryThreshold(symmetryThreshold);
    	// perform the calculations
    	System.out.println("Calculating power...");
    	long startTime = System.currentTimeMillis();
    	List<Power> results = null;
    	try
    	{
    		results = calc.getPower(params);
    	}
    	catch (Exception e)
    	{
    		System.err.println("Error in calculating power: " + e.getMessage());
    	}
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
    	
    	// accumulate results and calculate maximum absolute deviation
    	for(int i = 0; i < results.size(); i++, sasResultsIndex++)
    	{
    		GLMMPower power = (GLMMPower) results.get(i);
    		GLMMPower simPower = (simResults != null ? (GLMMPower) simResults.get(i) : null);
    		GLMMPower sasPower = (sasPowers != null ? (GLMMPower) sasPowers.get(sasResultsIndex) : null);
    		double simDeviation = Double.NaN;	
    		double sasDeviation = Double.NaN;	
    		double sasLowerCIDeviation = Double.NaN;
    		double sasUpperCIDeviation = Double.NaN;
    		
			ConfidenceInterval sasCI = null;
			ConfidenceInterval calcCI = null;
			
    		if (simPower != null)
    		{
    			simDeviation = Math.abs(power.getActualPower() - simPower.getActualPower());
    			if (simDeviation > maxSimDeviation) maxSimDeviation = simDeviation;
    		}
    		if (sasPower != null)
    		{
    			sasDeviation = Math.abs(power.getActualPower() - sasPower.getActualPower());
    			if (sasDeviation > maxSasDeviation) maxSasDeviation = sasDeviation;
    			sasCI = sasPower.getConfidenceInterval();
    			calcCI = power.getConfidenceInterval();
    			if (sasCI != null && calcCI != null)
    			{
    				sasLowerCIDeviation = Math.abs(calcCI.getLowerLimit() - sasCI.getLowerLimit());
    				sasUpperCIDeviation = Math.abs(calcCI.getUpperLimit() - sasCI.getUpperLimit());
    				if (sasLowerCIDeviation > maxSaslowerCIDeviation) maxSaslowerCIDeviation = sasLowerCIDeviation;
    				if (sasUpperCIDeviation > maxSasUpperCIDeviation) maxSasUpperCIDeviation = sasUpperCIDeviation;
    			}
    		}
    		
    		checkerResults.add(new Result(power,
    			(sasPower != null ? sasPower.getActualPower(): Double.NaN),
    			(sasCI != null ? sasCI.getLowerLimit(): Double.NaN),
    			(sasCI != null ? sasCI.getUpperLimit(): Double.NaN),
    			sasDeviation,
    			sasLowerCIDeviation,
    			sasUpperCIDeviation,
    			(simPower != null ? simPower.getActualPower() : Double.NaN),
    			simDeviation));
    	}
    }
    
    public void reset()
    {
    	checkerResults.clear();
    	timer.reset();
    	maxSasDeviation = 0;
    	maxSimDeviation = 0;
    	maxSaslowerCIDeviation = 0;
    	maxSasUpperCIDeviation = 0;
    	sasResultsIndex = 0;
    }
    
    /**
     * Write the results to an HTML file
     */
    public void outputResults(String title, String filename)
    {
    	outputResults(title, filename, null);
    }
    
    /**
     * Write the results to an HTML file, optionally including a plot image
     * @param filename
     */
    public void outputResults(String title, String filename, String imageFilename)
    {
    	// output the results
    	StringBuffer buffer = new StringBuffer();
    	
    	buffer.append("<html><head></head><body><h1>" + title + "</h1>");
    	buffer.append("<h3>Timing Results</h3>");
    	buffer.append("<table border='1' cellpadding='5'>");
    	buffer.append("<tr><td>Calculation time</td><td>" + ((double) timer.calculationMilliseconds / 1000)
    			+ "</td></tr>");
    	buffer.append("<tr><td>Simulation time</td><td>" + ((double) timer.simulationMilliseconds / 1000)
    			+ "</td></tr></table><p></p>");
    	
    	buffer.append("<h3>Max Absolute Deviation Summary</h3>");
    	buffer.append("<table border='1' cellpadding='5'>");
    	buffer.append("<tr><td>Max deviation from SAS</td><td>"+LongNumber.format(maxSasDeviation)+"</td></tr>");
    	if (maxSaslowerCIDeviation >= 0 && maxSasUpperCIDeviation >= 0)
    	{
    		buffer.append("<tr><td>Max deviation from lower CI limit</td><td>"+
    				LongNumber.format(maxSaslowerCIDeviation)+"</td></tr>");
    		buffer.append("<tr><td>Max deviation from upper CI limit</td><td>"+
    				LongNumber.format(maxSasUpperCIDeviation)+"</td></tr>");
    	}
    	buffer.append("<tr><td>Max deviation from simulation</td><td>"+
    			LongNumber.format(maxSimDeviation)+"</td></tr>");
    	buffer.append("</table><p></p>");
    	
    	buffer.append("<h3>Full Results</h3>");
    	buffer.append("<table border='1' cellpadding='5'><tr><th>Calc Power</th>");
    	if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
    	{
    		ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
    		buffer.append("<th>Confidence Interval (CI)</th><th>&alpha;-lower</th><th>&alpha;-upper</th>");
    	}
    	buffer.append("<th>SAS Power (deviation)</th>");
    	if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
    	{
    		ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
    		buffer.append("<th>SAS CI {deviation}</th>");
    	}
    	buffer.append("<th>Sim Power (deviation)</th>");
    	buffer.append("<th>Test</th><th>Sigma Scale</th><th>Beta Scale</th><th>Total N</th>");
    	buffer.append("<th>Alpha</th><th>Power Method</th><th>Quantile</th></tr>");

    	for(Result result: checkerResults)
    	{
    		buffer.append("<tr><td>");
    		buffer.append(Number.format(result.calculatedPower.getActualPower()));
    		buffer.append("</td><td>");
    		ConfidenceInterval ci = result.calculatedPower.getConfidenceInterval();
    		if (ci != null)
    		{
    			buffer.append("(" + Number.format(ci.getLowerLimit()) + ", " + Number.format(ci.getUpperLimit()) + ")");
        		buffer.append("</td><td>");
        		buffer.append(ci.getAlphaLower());
        		buffer.append("</td><td>");
        		buffer.append(ci.getAlphaUpper());
        		buffer.append("</td><td>");
    		}
    		buffer.append(Number.format(result.sasPower));
    		if (result.sasDeviation > sasTolerance)
    			buffer.append(" <font color='red'>(" + Number.format(result.sasDeviation) + ")</font></td><td>");
    		else
    			buffer.append(" (" + Number.format(result.sasDeviation) + ")</td><td>");
    		if (ci != null)
    		{
    			buffer.append("(" + Number.format(result.sasCILower) + ", " + Number.format(result.sasCIUpper) + ")");
        		buffer.append("{"); 
        		if (result.sasCILowerDeviation > sasTolerance)
        			buffer.append("<font color='red'>" + Number.format(result.sasCILowerDeviation) + "</font>");
        		else
        			buffer.append(Number.format(result.sasCILowerDeviation));
        		buffer.append(", ");
        		if (result.sasCIUpperDeviation > sasTolerance)
        			buffer.append("<font color='red'>" + Number.format(result.sasCIUpperDeviation) + "</font>");
        		else
        			buffer.append(Number.format(result.sasCIUpperDeviation));
        		buffer.append("}</td><td>");
    		}
    		buffer.append(Number.format(result.simulatedPower));
    		if (result.simulationDeviation > simTolerance)
    			buffer.append(" <font color='red'>(" + Number.format(result.simulationDeviation) + ")</font></td><td>");
    		else
    			buffer.append(" (" + Number.format(result.simulationDeviation) + ")</td><td>");

    		buffer.append(result.calculatedPower.getTest() + "</td><td>" + 
    				Number.format(result.calculatedPower.getSigmaScale()) + "</td><td>" + 
    				Number.format(result.calculatedPower.getBetaScale()) + "</td><td>" + 
    				result.calculatedPower.getTotalSampleSize() + "</td><td>" + 
    				Number.format(result.calculatedPower.getAlpha()) + "</td><td>" + 
    				powerMethodToString(result.calculatedPower.getPowerMethod()) + "</td><td>" + 
    				result.calculatedPower.getQuantile() + "</td></tr>");
    	}
    	
    	buffer.append("</table><p>");
    	
    	if (imageFilename != null)
    	{
    		buffer.append("<p><img src='" + imageFilename + "' /></p>");
    	}
    	
    	buffer.append("</body></html>");

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
    	System.out.println("Calc Power (lower, upper)\tSAS Power (deviation)\tSim Power (deviation)\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");
    	
    	for(Result result: checkerResults)
    	{
    		System.out.println(Number.format(result.calculatedPower.getActualPower()) + "(" + 
    				(result.calculatedPower.getConfidenceInterval() != null ? 
    						Number.format(result.calculatedPower.getConfidenceInterval().getLowerLimit()) : "n/a") + ", " + 
    				(result.calculatedPower.getConfidenceInterval() != null ? 
    						Number.format(result.calculatedPower.getConfidenceInterval().getUpperLimit()) : "n/a") + ")\t" +  
    	    		Number.format(result.sasPower) + " (" +  Number.format(result.sasDeviation) + ")\t" +
    				Number.format(result.simulatedPower) + " (" + Number.format(result.simulationDeviation) + ")\t" +
    				result.calculatedPower.getTest() + "\t" + 
    				Number.format(result.calculatedPower.getSigmaScale()) + "\t" + 
    				Number.format(result.calculatedPower.getBetaScale()) + "\t" + 
    				result.calculatedPower.getTotalSampleSize() + "\t" + 
    				Number.format(result.calculatedPower.getAlpha()) + "\t" + 
    				result.calculatedPower.getPowerMethod() + "\t" + 
    				result.calculatedPower.getQuantile() + "\t");
    	}
    	
    	System.out.println("Max Deviation from SAS: " + LongNumber.format(maxSasDeviation));
    	System.out.println("Max Deviation from Simulation: " + LongNumber.format(maxSimDeviation));
    	if (maxSaslowerCIDeviation >= 0 || maxSasUpperCIDeviation >=0)
    	{
        	System.out.println("Max Deviation from SAS, Lower Confidence Limit: " + 
        			LongNumber.format(maxSaslowerCIDeviation));
        	System.out.println("Max Deviation from SAS, Upper Confidence Limit: " + 
        			LongNumber.format(maxSasUpperCIDeviation));
    	}
    }
    
    private String powerMethodToString(GLMMPowerParameters.PowerMethod method)
    {
    	String value = null;
    	switch (method)
    	{
    	case CONDITIONAL_POWER:
    		value = "cond";
    		break;
    	case UNCONDITIONAL_POWER:
    		value = "uncond";
    		break;
    	case QUANTILE_POWER:
    		value = "quantile";
    		break;    		
    	}
    	return value;
    }
    
    public boolean isSASDeviationBelowTolerance()
    {
    	return (maxSasDeviation <= sasTolerance);
    }
    
    public boolean isSimulationDeviationBelowTolerance()
    {
    	return (maxSimDeviation <= simTolerance);
    }
    
    public void setSASTolerance(double tolerance)
    {
    	sasTolerance = tolerance;
    }
    
    public void setSimulationTolerance(double tolerance)
    {
    	simTolerance = tolerance;
    }
    
    public void setPositivityThreshold(double threshold)
    {
    	this.positivityThreshold = threshold;
    }
    
    public void setSymmetryThreshold(double threshold)
    {
    	this.symmetryThreshold = threshold;
    }
}
