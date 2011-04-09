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

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;

import junit.framework.TestCase;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;

/**
 * Unit test for fixed multivariate design with comparison against
 * simulation and SAS output.
 * @author Sarah Kreidler
 *
 */
public class TestConditionalMultivariate extends TestCase
{
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0,0.5,1,1.5,2};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final int[] SAMPLE_SIZE_LIST = {5,10};

	private static final String DATA_FILE =  "data" + File.separator + "TestConditionalMultivariate.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "FixedMultivariateOutput.html";
	private static final String TITLE = "Power results for fixed multivariate";
	
	private PowerChecker checker;
	
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
     * Test case for a multivariate GLM with only fixed predictors
     */
    public void testValidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();

        System.out.println(TITLE);
        checker.checkPower(params);
		checker.outputResults();
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
		assertTrue(checker.isSimulationDeviationBelowTolerance());
		checker.reset();
    }
    
    /**
     * Builds matrices for a multivariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildValidMultivariateFixedInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
     
        // add tests
        for(Test test: Test.values()) 
        {
            params.addTest(test);
        }

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int Q = 4;
        // create design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
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

}
