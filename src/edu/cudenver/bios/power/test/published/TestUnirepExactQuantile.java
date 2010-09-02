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
package edu.cudenver.bios.power.test.published;

import org.apache.commons.math.linear.Array2DRowRealMatrix;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RandomColumnMetaData;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker;
import junit.framework.TestCase;

/**
 * Test case for exact quantile power for the UNIREP.  Values should match
 * exact median power values from Table II in Glueck & Muller 2003.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestUnirepExactQuantile extends TestCase 
{
	private static final double MEAN = 9.75;
	private static final double VARIANCE = 2.0;
	private static final double[] ALPHA_LIST = {0.05};    
	private static final double[] SIGMA_SCALE_LIST = {1};	

	/**
	 * Compare the calculated UNIREP exact quantile powers against simulation
	 */
	public void testPower()
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

		System.out.println("Testing Multivariate, Random, Quantile");
		int mismatches = PowerChecker.checkPower(params5, false, false);
		mismatches += PowerChecker.checkPower(params25, false, false);
		mismatches += PowerChecker.checkPower(params50, false, false);

		assertEquals(0, mismatches);
	}

	/**
	 * Builds matrices for a multivariate GLM with a baseline covariate
	 * Note, this matrix set matches the values produced in Table II from Glueck&Muller
	 */
	private GLMMPowerParameters buildValidMultivariateRandomInputs(double[] betaScaleList, int repn)
	{
		GLMMPowerParameters params = new GLMMPowerParameters();
		params.setNonCentralityCDFExact(true);

		// add quantile power methods and median quantile
		params.addPowerMethod(PowerMethod.QUANTILE_POWER);
		params.addQuantile(0.5);

		// add UNIREP as the statistical test
		params.addTest(GLMMPowerParameters.Test.UNIREP);
		params.addTest(GLMMPowerParameters.Test.UNIREP_BOX);
		params.addTest(GLMMPowerParameters.Test.UNIREP_GEISSER_GREENHOUSE);
		params.addTest(GLMMPowerParameters.Test.UNIREP_HUYNH_FELDT);

		// add alpha values
		for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

		int P = 3;
		int Q = 3;
		// create design matrix
		double[][] essFixedData = {{1,0,0},{0,1,0},{0,0,1}};
		RowMetaData[] rowMd = {
				new RowMetaData(1), 
				new RowMetaData(1), 
				new RowMetaData(1)
		};
		double[][] essRandomData = {{1},{1},{1}};
		RandomColumnMetaData[] randColMd = {new RandomColumnMetaData(MEAN, VARIANCE)};
		DesignEssenceMatrix essence = new DesignEssenceMatrix(essFixedData, rowMd, 
				essRandomData, randColMd);
		params.setDesignEssence(essence);
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
		double [][] between = {{1,-1,0}, {1,0,-1}};
		double[][] betweenRandom = {{0}, {0}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));

		// build within subject contrast
		double [][] within = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

		return params;     
	}
}
