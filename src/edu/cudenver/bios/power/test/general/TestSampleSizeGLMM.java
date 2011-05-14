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

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

/**
 * Non-automated test for GLMM sample size
 * @author Sarah Kreidler
 *
 */
public class TestSampleSizeGLMM extends TestCase
{
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0,1,2};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final double[] POWER_LIST = {0.01,0.7,0.8,0.9};
    private DecimalFormat Number = new DecimalFormat("#0.000");
    
    public void testValidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
                
        List<Power> results = calc.getSampleSize(params);
        checkSampleSize(false, results);
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
        	List<Power> results = calc.getSampleSize(params);
        	checkSampleSize(false, results);
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
        
        List<Power> results = calc.getSampleSize(params);
        checkSampleSize(false, results);
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
                List<Power> results = calc.getSampleSize(params);
                checkSampleSize(false, results);
                fail();
        }
        catch (Exception e)
        {
            System.out.println("Correctly caught exception: " + e.getMessage());
        }
    }

    public void testValidMultivariateRandom()
    {
		// build the inputs
		double[] beta5 = {
				0.4997025,
				0.8075886,
				1.097641};
		GLMMPowerParameters params5 = buildValidMultivariateRandomInputs(beta5, 5);
		double[] beta25 = {
				0.1651525,
				0.2623301, 
				0.3508015
		};
		GLMMPowerParameters params25 = buildValidMultivariateRandomInputs(beta25, 25);
		double[] beta50 = {
				0.1141548,
				0.1812892,
				0.2423835
		};
		GLMMPowerParameters params50 = buildValidMultivariateRandomInputs(beta50, 50);

		GLMMPowerCalculator calc = new GLMMPowerCalculator();
		List<Power> sampleSizeList = calc.getSampleSize(params50); 
        checkSampleSize(true, sampleSizeList);
    }

    public void testInvalidMultivariateRandom()
    {

    }

    private void checkSampleSize(boolean isGlmmFG, List<Power> results)
    {
        int invalidSampleSize = 0;
        for(Power p: results)
        {
        	if (p.getActualPower() < p.getNominalPower()) invalidSampleSize++;
        }
        outputResults(isGlmmFG, results);
        assertEquals(invalidSampleSize, 0);
        
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
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(2));
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
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
//        for(Test test: Test.values()) 
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
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
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
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
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
     */
    private GLMMPowerParameters buildValidMultivariateRandomInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add power methods
        //for(PowerMethod method: PowerMethod.values()) params.addPowerMethod(method);
        params.addPowerMethod(PowerMethod.CONDITIONAL_POWER);
        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
        params.addQuantile(0.25);
        params.addQuantile(0.5);
        params.addQuantile(0.75);
        
        // add tests - only HL andUNIREP value for random case
        params.addTest(Test.HOTELLING_LAWLEY_TRACE);
        params.addTest(Test.UNIREP);
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int P = 3;
        int Q = 3;
        // create design matrix
        params.setDesignEssence(org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(Q));
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
        // build sigma G matrix
        double[][] sigmaG = {{VARIANCE}};
        params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

        // build sigma Y matrix
        double rho = 0.4;
        double [][] sigmaY = {{1,0},{0,1}};
        params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

        // build sigma YG
        double rhoYG = 0.8;
        double [][] sigmaYG = {{0.9},{0}};
        params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0},{0,0},{0,0}};
        double [][] betaRandom = {{1},{1}};
        params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,0,0}};
        double[][] betweenRandom = {{1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));
        
        // build within subject contrast
        double [][] within = {{1,0},{0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;     
    }
    
	/**
	 * Builds matrices for a multivariate GLM with a baseline covariate
	 * Note, this matrix set matches the values produced in Table II from Glueck&Muller
	 */
	private GLMMPowerParameters buildValidMultivariateRandomInputs(double[] betaScaleList, int repn)
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
		//  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
		params.addSampleSize(repn);
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
		// add beta scale values
		for(double betaScale: betaScaleList) params.addBetaScale(betaScale);

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

