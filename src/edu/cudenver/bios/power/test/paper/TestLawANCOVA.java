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
package edu.cudenver.bios.power.test.paper;

import java.text.DecimalFormat;
import java.util.List;

import junit.framework.TestCase;

import org.apache.commons.math.linear.Array2DRowRealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

/**
 * Test case to reproduce the power results from the GLIMMPSE manuscript 
 * for the ANCOVA example based on 
 * 
 * Law A, Logan H, Baron RS (1994). Desire for Control, Felt Control, and Stress Inoculation
 * Training during Dental Treatment." Journal of Personality and Social Psychology, 67(5),
 * 926-936.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestLawANCOVA extends TestCase
{
    private static DecimalFormat Number = new DecimalFormat("#0.0000");

	/**
	 * Reproduce the list of power displayed in the figure associated
	 * with the Law ANCOVA example in the GLIMMPSE manuscript
	 */
	public void testANCOVA()
	{
		try
		{
	        // build the inputs
	        GLMMPowerParameters params = new GLMMPowerParameters();
	        	
	        // add alpha values
	        params.addAlpha(0.01);
	        
	        // build the design matrix - reference cell coded
	        double [][] design = {
	        		{1,0,0,0,0,0,0,0},
	        		{1,0,0,1,0,0,0,0},
	        		{1,0,1,0,0,0,0,0},
	        		{1,0,1,1,0,0,1,0},
	        		{1,1,0,0,0,0,0,0},
	        		{1,1,0,1,0,1,0,0},
	        		{1,1,1,0,1,0,0,0},
	        		{1,1,1,1,1,1,1,1}
	        };
	        params.setDesignEssence(new Array2DRowRealMatrix(design));
	        
	        // add per group sample size
	        params.addSampleSize(20);
	        
	        // build beta matrix
	        double [][] beta = {{0},{0},{0},{0},{0},{0},{0},{1}};
	        double [][] betaRand = {{1}};
	        params.setBeta(new FixedRandomMatrix(beta, betaRand, false));
	        // add beta scale values
	        params.addBetaScale(2);
	        params.addBetaScale(4);
	        params.addBetaScale(6);
	        params.addBetaScale(8);
	        
	        // build between subject contrast
	        double [][] between = {{0,0,0,0,0,0,0,1}};
	        double [][] betweenRand = {{0}};
	        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRand, true));
	        
	        // build within subject contrast
	        double [][] within = {{1}};
	        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));
	        	
	        // build theta null matrix
	        double [][] theta0 = {{0}};
	        params.setTheta(new Array2DRowRealMatrix(theta0));
	        
	        // build covariance of Y
	        double [][] sigmaY = {{66}}; 
	        params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));
	        
	        // build covariance of g
	        double [][] sigmaG = {{66}}; 
	        params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));
	        
	        // build covariance of Y with g.  assumes correlation of 0.4 between Y and g
	        // sigmaYG = 0.4 * sqrt(66*66) = 26
	        double [][] sigmaYg = {{26}}; 
	        params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYg));
	        
	        // add sigma scale factors
	        params.addSigmaScale(0.5);
	        params.addSigmaScale(1);
	        params.addSigmaScale(2);
	        
	        // add tests
	        params.addTest(Test.HOTELLING_LAWLEY_TRACE);
	        
	        // select median quantile power
	        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
	        params.addQuantile(0.5);

	        // get the power results
	        GLMMPowerCalculator calc = new GLMMPowerCalculator();
	        List<Power> results = calc.getPower(params);
	        
	        // output the sample size results
			writeResults(results);
			assertTrue(true);			
		}
		catch (Exception e)
		{
			System.err.print(e.getMessage());
			fail();
		}
	}
	
	/**
	 * 
	 * @param results list of power results
	 */
	private void writeResults(List<Power> results)
	{
		System.out.println("Test, Actual Power, Total Sample Size," +
				"Beta Scale, Sigma Scale, Alpha, Nominal Power, Power Method, Quantile");

		for(Power power: results)
		{
			GLMMPower glmmPower = (GLMMPower) power;
			System.out.println(testToString(glmmPower.getTest()) + ", " + 
					Number.format(glmmPower.getActualPower()) + ", " + 
					glmmPower.getTotalSampleSize() + ", " + 
					glmmPower.getBetaScale() + ", " + 
					glmmPower.getSigmaScale() + ", " + 
					Number.format(glmmPower.getNominalPower()) + ", " + 
					powerMethodToString(glmmPower.getPowerMethod()) + ", " +
					glmmPower.getQuantile());
		}
	}
	
	/**
	 * Pretty print test name
	 */
	private String testToString(Test test)
	{
		switch (test)
		{
			case HOTELLING_LAWLEY_TRACE:
				return "hlt";
			case WILKS_LAMBDA:
				return "wl";
			case PILLAI_BARTLETT_TRACE:
				return "pbt";
			case UNIREP:
				return "unirep";
			case UNIREP_GEISSER_GREENHOUSE:
				return "unirepGG";
			case UNIREP_BOX:
				return "unirepBox";
			case UNIREP_HUYNH_FELDT:
				return "unirepHF";				
		}
		return "";
	}
	
	/**
	 * Pretty print power method
	 */
	private String powerMethodToString(PowerMethod method)
	{
		switch (method)
		{
		case CONDITIONAL_POWER:
			return "conditional";
		case UNCONDITIONAL_POWER:
			return "unconditional";
		case QUANTILE_POWER:
			return "quantile";
		}
		return "";
	}
}	
