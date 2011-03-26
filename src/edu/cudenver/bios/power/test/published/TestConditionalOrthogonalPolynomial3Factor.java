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
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import junit.framework.TestCase;

/**
 * Unit test for polynomial trends in U matrix
 * with comparison against simulation and SAS output.
 * Based on example 7 from POWERLIB (Johnson et al., 2009)
 * @author Sarah Kreidler
 *
 */
public class TestConditionalOrthogonalPolynomial3Factor extends TestCase
{
	private static final String DATA_FILE =  "sas" + File.separator + "data" + File.separator + "TestConditionalOrthogonalPolynomial3Factor.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestConditionalOrthogonalPolynomial3Factor.html";
	private static final String TITLE = "Power results for 3-way orthogonal polynomial contrasts";
	private PowerChecker checker;

	// groups for factors A,B, and C
	double[] factorA = {1,2,3};
	String nameA = "A";
	double[] factorB = {1,2,3};
	String nameB = "B";
	double[] factorC = {1,2,3};
	String nameC = "C";

	// times for factors D, E, and F
	double[] factorD = {1,2,3};
	String nameD = "D";
	double[] factorE = {1,2,3};
	String nameE = "E";
	double[] factorF = {1,2,3};
	String nameF = "F";

	public void setUp()
	{
		try
		{
			checker = new PowerChecker(DATA_FILE, false);
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			fail();
		}    
	}

	/**
	 * Test GLMM(F) with 1,2, and 3 Factor polynomial 
	 * contrasts in C and U matrices
	 */
	public void testOneToThreeFactorPolynomialContrasts()
	{
		System.out.println(TITLE);

		// add all tests
		//for(Test test: Test.values()) 
		Test[] testList = {Test.UNIREP};
		for(Test test: testList)
		{
			// set up the matrices
			GLMMPowerParameters params = buildInputsWithoutContrasts();
			params.addTest(test);
			// calculate all of the 3 factor polynomial contrasts
			OrthogonalPolynomialContrastCollection withinSubjectContrasts = 
				OrthogonalPolynomials.withinSubjectContrast(factorD, nameD, factorE, nameE, factorF, nameF);
			OrthogonalPolynomialContrastCollection betweenSubjectContrasts = 
				OrthogonalPolynomials.withinSubjectContrast(factorA, nameA, factorB, nameB, factorC, nameC);

			// run power for 1 factor contrasts		
			RealMatrix U = withinSubjectContrasts.getMainEffectContrast(nameD);
			RealMatrix C = betweenSubjectContrasts.getMainEffectContrast(nameA).transpose();
			RealMatrix thetaNull = MatrixUtils.createRealMatrixWithFilledValue(C.getRowDimension(),U.getColumnDimension(), 0);
			params.setWithinSubjectContrast(U);
			params.setBetweenSubjectContrast(new FixedRandomMatrix(C.getData(), null, true));
			params.setTheta(thetaNull);
			checker.checkPower(params);

//			// run power for 2 factor contrasts
//			U = withinSubjectContrasts.getTwoFactorInteractionContrast(nameD, nameE);
//			C = betweenSubjectContrasts.getTwoFactorInteractionContrast(nameA, nameB).transpose();
//			thetaNull = MatrixUtils.createRealMatrixWithFilledValue(C.getRowDimension(),U.getColumnDimension(), 0);
//			params.setWithinSubjectContrast(U);
//			params.setBetweenSubjectContrast(new FixedRandomMatrix(C.getData(), null, true));
//			params.setTheta(thetaNull);
//			checker.checkPower(params);
//
//			// run power for 3 factor contrasts
//			U = withinSubjectContrasts.getThreeFactorInteractionContrast(nameD, nameE, nameF);
//			C = betweenSubjectContrasts.getThreeFactorInteractionContrast(nameA, nameB, nameC).transpose();
//			thetaNull = MatrixUtils.createRealMatrixWithFilledValue(C.getRowDimension(),U.getColumnDimension(), 0);
//			params.setWithinSubjectContrast(U);
//			params.setBetweenSubjectContrast(new FixedRandomMatrix(C.getData(), null, true));
//			params.setTheta(thetaNull);
//			checker.checkPower(params);
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
	private GLMMPowerParameters buildInputsWithoutContrasts()
	{    	
		int Q = factorD.length * factorE.length * factorF.length;
		int P = factorA.length * factorB.length * factorC.length;
		// build the inputs        
		GLMMPowerParameters params = new GLMMPowerParameters();
		
		// add alpha values - bonferroni corrected for 6 comparisons
		params.addAlpha(0.05);

		// build the design matrix 
		params.setDesignEssence(org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(Q));

		// build beta matrix
		RealMatrix beta = MatrixUtils.createRealMatrixWithFilledValue(Q, P, 0);
		beta.setEntry(0, 0, 1);
		params.setBeta(new FixedRandomMatrix(beta.getData(), null, false));
		// add beta scale values
		params.addBetaScale(9);
		params.addBetaScale(18);
		params.addBetaScale(27);

		// build theta null matrix
		double [][] theta0 = {{0,0,0,0}};
		params.setTheta(new Array2DRowRealMatrix(theta0));

		// build sigma matrix
		RealMatrix sigma =org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(P);
		for(int i = 0; i < P; i++) sigma.setEntry(i, i, i+1);
		params.setSigmaError(sigma);
		// add sigma scale values
		params.addSigmaScale(1);

		// add sample size multipliers
		for(int perGroupN = 2; perGroupN <= 12; perGroupN += 2)
			params.addSampleSize(perGroupN);
		//params.addSampleSize(2);
		return params;
	}
}
