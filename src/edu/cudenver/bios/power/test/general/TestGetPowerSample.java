package edu.cudenver.bios.power.test.general;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RandomColumnMetaData;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

import junit.framework.TestCase;

public class TestGetPowerSample extends TestCase
{
	private static final double MEAN = 0;
	private static final double VARIANCE = 2.0;
	
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestGetPowerSample.dat";
	private static final String TITLE = "Power sample for a 2 sample t-test with a baseline covariate";
	
	public void testGetPowerSample()
	{
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();

        System.out.println(TITLE);
        // get the power sample
    	GLMMPowerCalculator calc = new GLMMPowerCalculator();
		List<double[]> powerSamples = calc.getSimulatedPowerSample(params, 1000);
		
		FileWriter writer = null;
		BufferedWriter out = null;
    	try
    	{
    		double[] powerValues = powerSamples.get(0);
    		System.out.println("Writing Output file with " + powerValues.length + " power values");
    		writer = new FileWriter(OUTPUT_FILE);
    		out = new BufferedWriter(writer);
    		int i = 0;
    		for(double value: powerValues) 
    		{
    			out.write(value + "\n");
    		}
    	} 
    	catch (IOException e) 
    	{
    		System.err.println("Failed to write output file");
    	}
    	finally
    	{
    			try { if (out != null) out.close();} catch (Exception e) {}
    			try {if (writer != null) writer.close(); } catch (Exception e) {}
    	}    	
	}
	
	/**
	 * Builds matrices for a multivariate GLM with a baseline covariate
	 * Note, this matrix set matches the values produced in Table II from Glueck&Muller
	 */
	private GLMMPowerParameters buildValidMultivariateRandomInputs()
	{
		GLMMPowerParameters params = new GLMMPowerParameters();

		// add quantile power method
		params.addPowerMethod(PowerMethod.QUANTILE_POWER);
		params.addQuantile(0.5);
		
		// add UNIREP as the statistical test
		params.addTest(Test.UNIREP);

		// add alpha values
		params.addAlpha(0.05);

		int P = 3;
		int Q = 3;
		// create design matrix
		params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
		// add sample size multipliers
		//  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
		params.addSampleSize(10);
		// build sigma G matrix
		double[][] sigmaG = {{VARIANCE}};
		params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

		// build sigma Y matrix
		double rho = 0.4;
		double [][] sigmaY = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

		// build sigma YG
		double [][] sigmaYG = {{0.5},{0.5}, {0.5}, {0}};
		params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

		// add sigma scale values
		params.addSigmaScale(1);

		// build beta matrix
		double [][] beta = {{1,0,0,0},{0,2,0,0},{0,0,0,0}};
		double [][] betaRandom = {{1,1,1,1}};
		params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));
		// add beta scale values
		params.addBetaScale(0.5);

		// build theta null matrix
		double [][] theta0 = {{0,0,0,0},{0,0,0,0}};
		params.setTheta(new Array2DRowRealMatrix(theta0));

		// build between subject contrast
		double [][] between = {{1,-1,0}, {1,0,-1}};
		double[][] betweenRandom = {{0}, {0}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));

		// build within subject contrast
		double [][] within = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

		return params;     
	}
	
}
