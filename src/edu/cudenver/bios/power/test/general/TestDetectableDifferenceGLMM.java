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

import junit.framework.TestCase;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

/**
 * Non-automated test cases for detectable difference feature of 
 * GLMMPowerCalculator.
 * @author Sarah Kreidler
 *
 */
public class TestDetectableDifferenceGLMM extends TestCase
{
	private static final double TOLERANCE = 1.0E-5;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final int[] SAMPLE_SIZE_LIST = {15,25,50};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final double[] POWER_LIST = {0.7,0.8,0.9};
    private DecimalFormat Number = new DecimalFormat("#0.000");

    
    public void testValidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
                
        List<Power> results = calc.getDetectableDifference(params);
        checkDetectableDifference(false, results);
    }

    public void testInvalidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);

        try
        {
        	List<Power> results = calc.getDetectableDifference(params);
            checkDetectableDifference(false, results);
            fail();
        }
        catch (Exception e)
        {
            System.out.println("Correctly caught exception: " + e.getMessage());
        }
    }

    public void testValidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        
        List<Power> results = calc.getDetectableDifference(params);
        checkDetectableDifference(false, results);

    }

    public void testInvalidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);

        try
        {
        	List<Power> results = calc.getDetectableDifference(params);
            checkDetectableDifference(false, results);
            fail();
        }
        catch (Exception e)
        {
            System.out.println("Exception: " + e.getMessage());
        }
    }

    public void testValidMultivariateRandom()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        
        List<Power> results = calc.getDetectableDifference(params);
        checkDetectableDifference(true, results);
    }

    
    private void checkDetectableDifference(boolean isGlmmFG, List<Power> results)
    {
        int invalidBetaScale = 0;
        for(Power p: results)
        {
        	if (Math.abs(p.getActualPower() - p.getNominalPower()) > TOLERANCE) invalidBetaScale++;
        }
        outputResults(isGlmmFG, results);
        assertEquals(invalidBetaScale, 0);
        
    }

    private void outputResults(boolean isGlmmFG, List<Power> resultList)
    {
    	System.out.println("Calc Power (lower, upper)\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");
    	
    	for(Power power: resultList)
    	{
    		GLMMPower result = (GLMMPower) power;
    		System.out.println(Number.format(result.getActualPower()) + "(" + 
    				(result.getConfidenceInterval() != null ? 
    						Number.format(result.getConfidenceInterval().getLowerLimit()) : "n/a") + ", " + 
    				(result.getConfidenceInterval() != null ? 
    						Number.format(result.getConfidenceInterval().getUpperLimit()) : "n/a") + ")\t" +  
    				result.getTest() + "\t" + 
    				Number.format(result.getSigmaScale()) + "\t" + 
    				Number.format(result.getBetaScale()) + "\t" + 
    				result.getTotalSampleSize() + "\t" + 
    				Number.format(result.getAlpha()) + "\t" + 
    				result.getPowerMethod() + "\t" + 
    				result.getQuantile() + "\t");
    	}
    }
    
    /********** helper functions to create the matrices ***********/
    
    /**
     * Builds matrices for a univariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildValidUnivariateInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
       
        // add tests
        for(Test test: Test.values()) 
        {
            params.addTest(test);
        }
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build design matrix
        double[][] essenceData = {{1,0},{0,1}};
        params.setDesignEssence(new Array2DRowRealMatrix(essenceData));
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        return params;     
    }   

    /**
     * Builds matrices for a multivariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildValidMultivariateFixedInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
     
        // add tests
//        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
//        {
//            params.addTest(test);
//        }
        params.addTest(Test.WILKS_LAMBDA);
        params.addTest(Test.PILLAI_BARTLETT_TRACE);
        params.addTest(Test.HOTELLING_LAWLEY_TRACE);

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int Q = 4;
        // create design matrix
        RealMatrix essenceMatrix = MatrixUtils.createRealIdentityMatrix(Q);
        params.setDesignEssence(essenceMatrix);
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build sigma matrix
        double rho = 0.4;
        double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}}; // compound symmetry
        // double [][] sigma = {{1,0.2,0.3},{0.2,1,0.2},{0.3,0.2,1}}; // toeplitz
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null,  false));
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        // build within subject contrast
        double [][] within = {{1,1},{-1,0},{0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        //RealMatrix U = params.getWithinSubjectContrast();
        //RealMatrix upu = U.multiply(U.transpose());
        
        return params;     
    }   
    
	/**
	 * Builds matrices for a multivariate GLM with a baseline covariate
	 * Note, this matrix set matches the values produced in Table II from Glueck&Muller
	 */
	private GLMMPowerParameters buildValidMultivariateRandomInputs()
	{
		GLMMPowerParameters params = new GLMMPowerParameters();

        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
		// add quantile power method
		params.addPowerMethod(PowerMethod.QUANTILE_POWER);
		params.addQuantile(0.5);
		
		// add HLT as the statistical test
		params.addTest(Test.HOTELLING_LAWLEY_TRACE);

		// add alpha values
		for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

		int P = 3;
		int Q = 3;
		// create design matrix
		params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
		// add sample size multipliers
		for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);

		// build sigma G matrix
		double[][] sigmaG = {{1}};
		params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

		// build sigma Y matrix
		double rho = 0.4;
		double [][] sigmaY = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

		// build sigma YG
		double [][] sigmaYG = {{0.5},{0.5}, {0.5}, {0}};
		params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

		// add sigma scale values
		for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);

		// build beta matrix
		double [][] beta = {{1,0,0,0},{0,2,0,0},{0,0,0,0}};
		double [][] betaRandom = {{1,1,1,1}};
		params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));

		// build theta null matrix
		double [][] theta0 = {{0,0,0,0},{0,0,0,0}};
		params.setTheta(new Array2DRowRealMatrix(theta0));

		// build between subject contrast
		double [][] between = {{-1,1,0}, {-1,0,1}};
		double[][] betweenRandom = {{0}, {0}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));

		// build within subject contrast
		double [][] within = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

		return params;     
	}
    
}

