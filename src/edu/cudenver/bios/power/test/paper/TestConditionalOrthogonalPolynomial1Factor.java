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

import java.io.File;
import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
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
public class TestConditionalOrthogonalPolynomial1Factor extends TestCase
{
	private static final String DATA_FILE =  "data" + File.separator + "TestConditionalOrthogonalPolynomial1Factor.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestConditionalOrthogonalPolynomial1Factor.html";
	private static final String TITLE = "GLMM(F) Example 7. Power for a time by treatment interaction using orthogonal polynomial contrast for time";
	private PowerChecker checker;
	private boolean verbose = false;
	
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
     * Test GLMM(F) with polynomial contrasts in U matrix
     */
    public void testPolynomial1Factor()
    {
        // build the inputs        
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests
        for(Test test: Test.values()) 
        {
            params.addTest(test);
        }
        
        // add alpha values - bonferroni corrected for 6 comparisons
        params.addAlpha(0.05);

        // build beta matrix
        double [][] beta = {{0,0,0,0,1},{1,0,0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        params.addBetaScale(1);
        
        // build theta null matrix
        double [][] theta0 = {{0,0,0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double rho = 0.375;
        double var = 1.5;
        double [][] sigma = {
        		{var,rho,rho,rho,rho},
        		{rho,var,rho,rho,rho},
        		{rho,rho,var,rho,rho},
        		{rho,rho,rho,var,rho},
        		{rho,rho,rho,rho,var}
        };
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        params.addSigmaScale(1);
        
        // build design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(2));
        // add sample size multipliers
        params.addSampleSize(10);
        params.addSampleSize(20);
        params.addSampleSize(40);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));
        
        // build within subject contrast
        double[] times ={2, 4 ,6, 8, 10};
        String name = "times";
        RealMatrix U = OrthogonalPolynomials.withinSubjectContrast(times, name).getMainEffectContrast(name);
        if (verbose) printMatrix("U Matrix", U);
        params.setWithinSubjectContrast(U);
        
        checker.checkPower(params);
		checker.outputResults(TITLE);
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
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
}
