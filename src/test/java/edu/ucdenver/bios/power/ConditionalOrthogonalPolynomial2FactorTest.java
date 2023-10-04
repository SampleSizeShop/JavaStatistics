/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2015 Regents of the University of Colorado.
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
package edu.ucdenver.bios.power;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.OrthogonalPolynomialContrastCollection;
import edu.cudenver.bios.matrix.OrthogonalPolynomials;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateEpsilonApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import edu.cudenver.bios.utils.Factor;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * Unit test for polynomial trends in U matrix
 * with comparison against simulation and SAS output.
 * Based on example 9 from POWERLIB (Johnson et al., 2009)
 * @author Sarah Kreidler
 *
 */
public class ConditionalOrthogonalPolynomial2FactorTest {

	private static final String MB_DATA_FILE =  "TestConditionalOrthogonalPolynomial2FactorMB.xml";
	private static final String MEST_DATA_FILE =  "TestConditionalOrthogonalPolynomial2FactorMEST.xml";

	private static final String TITLE_MB = "GLMM(F) Example 9 MB: Power for a multivariate model" +
			" with two within subject factors, using the Muller and Barton approximation";
	private static final String TITLE_MEST = "GLMM(F) Example 9 MEST: Power for a multivariate " +
			"model with two within subject factors, using the Muller, Edwards, Simpson, Taylor approximation";
	private static final double TOLERANCE = 0.000001;

	private boolean verbose = false;
	// groups for factors A,B, and C
	double[] dataA = {1,2,4};
	String nameA = "A";
	Factor factorA = new Factor(nameA, dataA);
	double[] dataB = {1,3,5};
	String nameB = "B";
	Factor factorB = new Factor(nameB, dataB);
	
	// sigma* matrices
	double[] varE = {.47960, .01000, .01000, .01000}; // epsilon ~ 0.28
	RealMatrix sigStarE = 
		org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(varE);
	double[] varF = {.34555, .06123, .05561, .04721}; // epsilon ~ 0.5
	RealMatrix sigStarF = 
		org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(varF);
	double[] varG =  {.23555, .17123, .05561, .04721}; // epsilon ~ .72 
	RealMatrix sigStarG = 
		org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(varG);
	double[] varH = {.12740, .12740, .12740, .12740}; // epsilon = 1 
	RealMatrix sigStarH = 
		org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(varH);
	RealMatrix[] sigmaStars = {sigStarE, sigStarF, sigStarG, sigStarH};

	/**
	 * Test GLMM(F) with 2 Factor polynomial U matrix using the
	 * Muller & Barton (1989) approximation
	 */
	@Test
	public void testTwoFactorContrastMullerBarton()
	{
		List<GLMMPower> sasPowers = Utils.readSasPowers(MB_DATA_FILE);
		PowerChecker checker = new PowerChecker(sasPowers, true);
		checker.setSymmetryThreshold(1.0E-14);

		// set up the matrices and approximation method
		GLMMPowerParameters params = buildSharedInputs();
		// calculate power using Muller & Barton approximations
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP,
				UnivariateCdfApproximation.MULLER_BARTON_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_BOX,
				UnivariateCdfApproximation.MULLER_BARTON_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE,
				UnivariateCdfApproximation.MULLER_BARTON_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_HUYNH_FELDT,
				UnivariateCdfApproximation.MULLER_BARTON_APPROX);
		// set epsilon method to Muller Barton
		params.setUnivariateEpsilonMethod(GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE,
				UnivariateEpsilonApproximation.MULLER_BARTON_APPROX);
		params.setUnivariateEpsilonMethod(GLMMTestFactory.Test.UNIREP_HUYNH_FELDT,
				UnivariateEpsilonApproximation.MULLER_BARTON_APPROX);

		checkPower(checker, TITLE_MB, params);
	}
	
	/**
	 * Test GLMM(F) with 2 Factor polynomial U matrix using the
	 * Muller & Barton (1989) approximation
	 */
	@Test
	public void testTwoFactorContrastMullerEdwardsSimpsonTaylor() {
		List<GLMMPower> sasPowers = Utils.readSasPowers(MEST_DATA_FILE);
		PowerChecker checker = new PowerChecker(sasPowers, true);
		checker.setSymmetryThreshold(1.0E-14);
		
		// set up the matrices and approximation method
		GLMMPowerParameters params = buildSharedInputs();
		// calculate power using Muller, Edwards, Taylor approximations
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP,
				UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_BOX,
				UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE,
				UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
		params.setUnivariateCdfMethod(GLMMTestFactory.Test.UNIREP_HUYNH_FELDT,
				UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
		// set epsilon method to Muller, Edwards, Taylor approximation
		params.setUnivariateEpsilonMethod(GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE,
				UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX);
		params.setUnivariateEpsilonMethod(GLMMTestFactory.Test.UNIREP_HUYNH_FELDT,
				UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX);

		checkPower(checker, TITLE_MEST, params);
	}

	/**
	 * 
	 * @param checker
	 * @param params
	 */
	private void checkPower(PowerChecker checker, String title, GLMMPowerParameters params)
	{
		/* 
		 * get orthogonal contrasts for within subject factors
		 * Log base 2 spacing Clip (2,4,16) and Region(2,8,32) 
		 */
        ArrayList<Factor> factorList = new ArrayList<Factor>();
        factorList.add(factorA);
        factorList.add(factorB);
		OrthogonalPolynomialContrastCollection collection = 
			OrthogonalPolynomials.withinSubjectContrast(factorList);
		params.setWithinSubjectContrast(collection.getInteractionContrast(factorList).getContrastMatrix());
		
		// theta critical matrix used to back-tranform beta from U
		//	 THETA = {.25}#{.5 1 -1 .5}; * =Theta(cr) from 1st sentence *after* 
		//  equation 7, Coffey and Muller (2003); 
		double [][] thetaData = {{0.5,1,-1,0.5}};
		RealMatrix thetaCr = (new Array2DRowRealMatrix(thetaData)).scalarMultiply(0.25);
		
		// loop over the various sigma matrices
		//Test[] testList = Test.values();
		GLMMTestFactory.Test[] testList = {
				GLMMTestFactory.Test.UNIREP,
				GLMMTestFactory.Test.UNIREP_BOX,
				GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE,
				GLMMTestFactory.Test.UNIREP_HUYNH_FELDT
		};
		for (GLMMTestFactory.Test test : testList) {
			params.clearTestList();
			params.addTest(test);
			
			for(RealMatrix sigStar: sigmaStars)
			{
				RealMatrix U = params.getWithinSubjectContrast();
				// 1st paragraph in section 2.4, Coffey and Muller 2003 *;
				RealMatrix sigmaTemp = U.multiply(sigStar).multiply(U.transpose());
				int dimension = sigmaTemp.getRowDimension();
				RealMatrix Uother = 
					MatrixUtils.getRealMatrixWithFilledValue(dimension, 1, 1/Math.sqrt(U.getColumnDimension()));
				Uother = MatrixUtils.getHorizontalAppend(Uother, collection.getMainEffectContrast(factorA).getContrastMatrix());
				Uother = MatrixUtils.getHorizontalAppend(Uother, collection.getMainEffectContrast(factorB).getContrastMatrix());
				double varianceMean = (double) sigStar.getTrace() / (double) sigStar.getColumnDimension();
				RealMatrix sigmaError = sigmaTemp.add(Uother.multiply(Uother.transpose()).scalarMultiply(varianceMean));
				
				if (verbose) printMatrix("Sigma Error",sigmaError);
				// 1st paragraph in section 2.4, Coffey and Muller 2003 *; 
				RealMatrix beta = thetaCr.multiply(U.transpose());       
				params.setSigmaError(sigmaError);
				params.setBeta(new FixedRandomMatrix(beta.getData(), null, false));

				checker.checkPower(params);
			}

		}

        // output the results
		ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
		reportBuilder.createValidationReportAsStdout(checker, title, true);

		assertTrue("results outside tolerance: " + TOLERANCE, checker.isSASDeviationBelowTolerance(TOLERANCE));
	}
	
	
	/**
	 * Write the matrix to std out
	 * @param m
	 */
	private void printMatrix(String title, RealMatrix m)
	{
		System.out.println(title);
		DecimalFormat Number = new DecimalFormat("#0.00000000000000000000");
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
		params.setDesignEssence(org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(1));

		// set between subject contrast
		double[][] cData = {{1}};
		params.setBetweenSubjectContrast(new FixedRandomMatrix(cData, null, true));
		
		// add beta scale values
		params.addBetaScale(1);

		// build theta null matrix
		params.setTheta(MatrixUtils.getRealMatrixWithFilledValue(1, 4, 0));

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
