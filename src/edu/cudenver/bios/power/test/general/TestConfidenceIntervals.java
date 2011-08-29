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

import java.io.File;

import junit.framework.TestCase;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;

/**
 * General unit test for confidence intervals.  Covers univariate and multivariate
 * with beta known and unknown 
 * 
 * @author Sarah Kreidler
 *
 */
public class TestConfidenceIntervals extends TestCase
{
	private static final String UNIVARIATE_DATA_FILE = "data" + 
	File.separator + "TestUnivariateConfidenceIntervals.xml";
	private static final String UNIVARIATE_OUTPUT_FILE = "text" + File.separator + "results" + 
	File.separator + "TestUnivariateConfidenceIntervals.html";
	private static final String MULTIVARIATE_BETA_KNOWN_DATA_FILE = "data" + 
	File.separator + "TestMultivariateConfidenceIntervalsBetaKnown.xml";
	private static final String MULTIVARIATE_BETA_KNOWN_OUTPUT_FILE = "text" + File.separator + "results" + 
	File.separator + "TestMultivariateConfidenceIntervalsBetaKnown.html";
	private static final String MULTIVARIATE_BETA_UNKNOWN_DATA_FILE = "data" + 
	File.separator + "TestMultivariateConfidenceIntervalsBetaUnknown.xml";
	private static final String MULTIVARIATE_BETA_UNKNOWN_OUTPUT_FILE = "text" + File.separator + "results" + 
	File.separator + "TestMultivariateConfidenceIntervalsBetaUnknown.html";

	/**
	 * Compare CI results between JavaStatistics 
	 * and POWERLIB
	 */
	public void testUnivariateWithConfidenceLimits()
	{
		try
		{
			PowerChecker checker;
			checker = new PowerChecker(UNIVARIATE_DATA_FILE, false);

			GLMMPowerParameters params = new GLMMPowerParameters();

			// add alpha values
			params.addAlpha(0.01);

			// build beta matrix
			double [][] beta = {{0},{1}};
			params.setBeta(new FixedRandomMatrix(beta, null, false));
			// beta scale of 1 and 0
			params.addBetaScale(0);
			params.addBetaScale(0.25);
			params.addBetaScale(0.5);
			params.addBetaScale(0.75);

			// build theta null matrix
			double [][] theta0 = {{0}};
			params.setTheta(new Array2DRowRealMatrix(theta0));

			// build sigma matrix
			double [][] sigma = {{0.068}};
			params.setSigmaError(new Array2DRowRealMatrix(sigma));
			// add sigma scale values
			params.addSigmaScale(1);

			// build design matrix
			params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(2));
			// add sample size multipliers
			params.addSampleSize(12);

			// build between subject contrast
			double [][] between = {{1,-1}};
			params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

			// parameters for confidence limits
			params.setSampleSizeForEstimates(24);
			params.setDesignMatrixRankForEstimates(2);

			// run the test for both types of CI's
			for(Test test: Test.values()) 
			{
				//			Test[] tmpTests = {Test.UNIREP};
				//			for(Test test: tmpTests)
				//			{
				params.clearTestList();
				params.addTest(test);

				for(ConfidenceIntervalType ciType: ConfidenceIntervalType.values())
				{
					if (ciType == ConfidenceIntervalType.NONE) continue;
					params.setConfidenceIntervalType(ciType);
					// 2 sided CI
					params.setAlphaLowerConfidenceLimit(0.025);
					params.setAlphaUpperConfidenceLimit(0.025);
					checker.checkPower(params);
					// 1 sided lower CI
					params.setAlphaLowerConfidenceLimit(0.05);
					params.setAlphaUpperConfidenceLimit(0);
					checker.checkPower(params);
					// 1 sided upper CI
					params.setAlphaLowerConfidenceLimit(0);
					params.setAlphaUpperConfidenceLimit(0.05);
					checker.checkPower(params);
				}
			}
			// output the results
			String title = "Univariate confidence intervals with beta known and unknown";
			checker.outputResults(title);
			checker.outputResults(title, UNIVARIATE_OUTPUT_FILE);
			assertTrue(checker.isSASDeviationBelowTolerance());
			checker.reset();
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			e.printStackTrace();
			fail();
		}
	}

	/**
	 * Compare multivariate CI results between JavaStatistics
	 * and POWERLIB when beta is assumed known
	 */
	public void testMultivariateWithConfidenceLimitsBetaKnown()
	{
		try
		{
			PowerChecker checker;
			checker = new PowerChecker(MULTIVARIATE_BETA_KNOWN_DATA_FILE, false);

			GLMMPowerParameters params = buildMultivariateInputs();

			// add all tests
			for(Test test: Test.values()) 
			{
				params.clearTestList();
				params.addTest(test);

				// calculate power, CI's assuming beta is known
				params.setConfidenceIntervalType(ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED);
				// 2 sided CI
				params.setAlphaLowerConfidenceLimit(0.025);
				params.setAlphaUpperConfidenceLimit(0.025);
				checker.checkPower(params);
				// 1 sided lower CI
				params.setAlphaLowerConfidenceLimit(0.05);
				params.setAlphaUpperConfidenceLimit(0);
				checker.checkPower(params);
				// 1 sided upper CI
				params.setAlphaLowerConfidenceLimit(0);
				params.setAlphaUpperConfidenceLimit(0.05);
				checker.checkPower(params);
			}
			
			// output the results
			String title = "Multivariate confidence intervals with beta known";
			checker.outputResults(title);
			checker.outputResults(title, MULTIVARIATE_BETA_KNOWN_OUTPUT_FILE);
			assertTrue(checker.isSASDeviationBelowTolerance());
			checker.reset();
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			e.printStackTrace();
			fail();
		}
	}
	
	
	/**
	 * Compare multivariate CI results between JavaStatistics
	 * and POWERLIB when beta is assumed known
	 */
	public void testMultivariateWithConfidenceLimitsBetaUnknown()
	{
		try
		{
			PowerChecker checker;
			checker = new PowerChecker(MULTIVARIATE_BETA_UNKNOWN_DATA_FILE, false);

			GLMMPowerParameters params = buildMultivariateInputs();

			// add all tests
			Test[] testList = {Test.WILKS_LAMBDA, Test.PILLAI_BARTLETT_TRACE, Test.HOTELLING_LAWLEY_TRACE};
			for(Test test: testList) 
			{
				params.clearTestList();
				params.addTest(test);

				// calculate power, CI's assuming beta is known
				params.setConfidenceIntervalType(ConfidenceIntervalType.BETA_SIGMA_ESTIMATED);
				// 2 sided CI
				params.setAlphaLowerConfidenceLimit(0.025);
				params.setAlphaUpperConfidenceLimit(0.025);
				checker.checkPower(params);
				// 1 sided lower CI
				params.setAlphaLowerConfidenceLimit(0.05);
				params.setAlphaUpperConfidenceLimit(0);
				checker.checkPower(params);
				// 1 sided upper CI
				params.setAlphaLowerConfidenceLimit(0);
				params.setAlphaUpperConfidenceLimit(0.05);
				checker.checkPower(params);
			}
			
			// output the results
			String title = "Multivariate confidence intervals with beta unknown";
			checker.outputResults(title);
			checker.outputResults(title, MULTIVARIATE_BETA_KNOWN_OUTPUT_FILE);
			assertTrue(checker.isSASDeviationBelowTolerance());
			checker.reset();
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			e.printStackTrace();
			fail();
		}
	}

	private GLMMPowerParameters buildMultivariateInputs()
	{
		GLMMPowerParameters params = new GLMMPowerParameters();

		// add alpha values
		params.addAlpha(0.01);

		// build design matrix
		params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(3));
		// add sample size multipliers
		params.addSampleSize(12);
		
		// build beta matrix
		double [][] beta = {{1,0,0},{0,0,0},{0,0,0}};
		params.setBeta(new FixedRandomMatrix(beta, null, false));
		// beta scales
		double[] betaScale = {0.75, 1, 2};
		for(double scale: betaScale) params.addBetaScale(scale);

		// build theta null matrix
		double [][] theta0 = {{0,0},{0,0}};
		params.setTheta(new Array2DRowRealMatrix(theta0));

		// build sigma matrix
		double rho = 0.4;
		double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}};
		params.setSigmaError(new Array2DRowRealMatrix(sigma));
		// add sigma scale values
		params.addSigmaScale(1);

		// build between subject contrast
		double [][] between = {{1,-1,0},{1,0,-1}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

		// build between subject contrast
		double [][] within = {{1,1}, {-1,0},{0,-1}};
		params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

		// parameters for confidence limits
		params.setSampleSizeForEstimates(24);
		params.setDesignMatrixRankForEstimates(2);
		
		return params;
	}
}
