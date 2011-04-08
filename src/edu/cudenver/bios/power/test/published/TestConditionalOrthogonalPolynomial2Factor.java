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

import java.io.File;
import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.OrthogonalPolynomialContrastCollection;
import edu.cudenver.bios.matrix.OrthogonalPolynomials;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateEpsilonApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import junit.framework.TestCase;

/**
 * Unit test for polynomial trends in U matrix
 * with comparison against simulation and SAS output.
 * Based on example 9 from POWERLIB (Johnson et al., 2009)
 * @author Sarah Kreidler
 *
 */
public class TestConditionalOrthogonalPolynomial2Factor extends TestCase
{
	private static final String DATA_FILE =  "sas" + File.separator + "data" + File.separator + "TestConditionalOrthogonalPolynomial2Factor.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestConditionalOrthogonalPolynomial2Factor.html";
	private static final String TITLE = "Power results for 2 within factor orthogonal polynomial contrasts";
	private PowerChecker checker;
	
	// groups for factors A,B, and C
	double[] factorA = {1,2,4};
	String nameA = "A";
	double[] factorB = {1,3,5};
	String nameB = "B";

	// sigma* matrices
	double[] varE = {.47960, .01000, .01000, .01000}; // epsilon ~ 0.28
	RealMatrix sigStarE = 
		org.apache.commons.math.linear.MatrixUtils.createRealDiagonalMatrix(varE);
	double[] varF = {.34555, .06123, .05561, .04721}; // epsilon ~ 0.5
	RealMatrix sigStarF = 
		org.apache.commons.math.linear.MatrixUtils.createRealDiagonalMatrix(varF);
	double[] varG =  {.23555, .17123, .05561, .04721}; // epsilon ~ .72 
	RealMatrix sigStarG = 
		org.apache.commons.math.linear.MatrixUtils.createRealDiagonalMatrix(varG);
	double[] varH = {.12740, .12740, .12740, .12740}; // epsilon = 1 
	RealMatrix sigStarH = 
		org.apache.commons.math.linear.MatrixUtils.createRealDiagonalMatrix(varH);
	RealMatrix[] sigmaStars = {sigStarE, sigStarF, sigStarG, sigStarH};
	
	public void setUp()
	{
		try
		{
			checker = new PowerChecker(DATA_FILE, true);
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			fail();
		}    
	}

	/**
	 * Test GLMM(F) with 2 Factor polynomial U matrix
	 */
	public void testTwoFactorPolynomialContrasts()
	{
		System.out.println(TITLE);
		
		// set up the matrices
		GLMMPowerParameters params = buildSharedInputs();
		
		// theta critical matrix used to back-tranform beta from U
		//	 THETA = {.25}#{.5 1 -1 .5}; * =Theta(cr) from 1st sentence *after* 
		//  equation 7, Coffey and Muller (2003); 
		double [][] thetaData = {{0.5,1,-1,0.5}};
		RealMatrix thetaCr = (new Array2DRowRealMatrix(thetaData)).scalarMultiply(0.25);
		
		// loop over the various sigma matrices
		//Test[] testList = Test.values();
		Test[] testList = {Test.UNIREP, Test.UNIREP_BOX, 
				Test.UNIREP_GEISSER_GREENHOUSE, Test.UNIREP_HUYNH_FELDT};
		for(Test test : testList)
		{
			params.clearTestList();
			params.addTest(test);
			
			for(RealMatrix sigStar: sigmaStars)
			{
				RealMatrix U = params.getWithinSubjectContrast();
				// 1st paragraph in section 2.4, Coffey and Muller 2003 *;
				RealMatrix sigmaError = U.multiply(sigStar).multiply(U.transpose());
				// 1st paragraph in section 2.4, Coffey and Muller 2003 *; 
				RealMatrix beta = thetaCr.multiply(U.transpose());       
				params.setSigmaError(sigmaError);
				params.setBeta(new FixedRandomMatrix(beta.getData(), null, false));
			
				// calculate power using Muller & Barton approximations
				params.setUnivariateCdfMethod(Test.UNIREP, 
						UnivariateCdfApproximation.MULLER_BARTON_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_BOX, 
						UnivariateCdfApproximation.MULLER_BARTON_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_GEISSER_GREENHOUSE, 
						UnivariateCdfApproximation.MULLER_BARTON_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_HUYNH_FELDT, 
						UnivariateCdfApproximation.MULLER_BARTON_APPROX);
				// set epsilon method to Muller Barton
				params.setUnivariateEpsilonMethod(Test.UNIREP_GEISSER_GREENHOUSE, 
						UnivariateEpsilonApproximation.MULLER_BARTON_APPROX);
				params.setUnivariateEpsilonMethod(Test.UNIREP_HUYNH_FELDT, 
						UnivariateEpsilonApproximation.MULLER_BARTON_APPROX);
				checker.checkPower(params);

				// calculate power using Muller, Edwards, Taylor approximations
				params.setUnivariateCdfMethod(Test.UNIREP, 
						UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_BOX, 
						UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_GEISSER_GREENHOUSE, 
						UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				params.setUnivariateCdfMethod(Test.UNIREP_HUYNH_FELDT, 
						UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				// set epsilon method to Muller, Edwards, Taylor approximation
				params.setUnivariateEpsilonMethod(Test.UNIREP_GEISSER_GREENHOUSE, 
						UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				params.setUnivariateEpsilonMethod(Test.UNIREP_HUYNH_FELDT, 
						UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
				checker.checkPower(params);
			}

		}
		// output and test the results
		checker.outputResults();
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
		assertTrue(checker.isSimulationDeviationBelowTolerance());
		checker.reset();
	}

	
	
	/**
	 * Write the matrix to std out
	 * @param m
	 */
	private void printMatrix(String title, RealMatrix m)
	{
		System.out.println(title);
		DecimalFormat Number = new DecimalFormat("#0.000");
		for(int row = 0; row < m.getRowDimension(); row++)
		{
			for(int col= 0; col < m.getColumnDimension(); col++)
			{
				System.out.print(Number.format(m.getEntry(row, col)) + "\t");
			}
			System.out.print("\n");
		}
	}

	/**
	 * Sets up the basic X, beta matrices for each of the tests.  The individual
	 * tests create different contrast matrices
	 * @return
	 */
	private GLMMPowerParameters buildSharedInputs()
	{    	

		// build the inputs        
		GLMMPowerParameters params = new GLMMPowerParameters();
		
		// add alpha values - bonferroni corrected for 6 comparisons
		params.addAlpha(0.04);

		// build the design matrix 
		params.setDesignEssence(org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(1));

		/* 
		 * get orthogonal contrasts for within subject factors
		 * Log base 2 spacing Clip (2,4,16) and Region(2,8,32) 
		 */
		OrthogonalPolynomialContrastCollection collection = 
			OrthogonalPolynomials.withinSubjectContrast(factorA, nameA, factorB, nameB);
		params.setWithinSubjectContrast(collection.getTwoFactorInteractionContrast(nameA, nameB));
		// set between subject contrast
		double[][] cData = {{1}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(cData, null, true));
		
		// add beta scale values
		params.addBetaScale(1);

		// build theta null matrix
		params.setTheta(MatrixUtils.createRealMatrixWithFilledValue(1, 4, 0));

		// add sigma scale values
		// gamma in Coffey and Muller (2003) *;
		params.addSigmaScale(0.5);
		params.addSigmaScale(1.0);
		params.addSigmaScale(2.0);
		
		// add sample size multipliers
		params.addSampleSize(20);

		return params;
	}
}