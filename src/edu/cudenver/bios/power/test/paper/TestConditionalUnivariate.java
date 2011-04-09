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

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import junit.framework.TestCase;

/**
 * Unit test for fixed univariate design with comparison against
 * simulation and SAS output.
 * @author Sarah Kreidler
 *
 */
public class TestConditionalUnivariate extends TestCase
{
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0,0.5,1,1.5,2};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final int[] SAMPLE_SIZE_LIST = {10};

	private static final String DATA_FILE =  "data" + File.separator + "TestConditionalUnivariate.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "FixedUnivariateOutput.html";
	private static final String TITLE = "Power results for fixed univariate";
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
     * Test valid inputs for a univariate linear model with only fixed predictors
     */
    public void testValidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        System.out.println(TITLE);
        checker.checkPower(params);
		checker.outputResults();
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
		assertTrue(checker.isSimulationDeviationBelowTolerance());
		checker.reset();
    }

    /**
     * Tests if the calculator throws an exception on invalid inputs
     */
    public void testInvalidInputs()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);
        
        try
        {
        	calc.getPower(params);
            fail();
        }
        catch (Exception e)
        {
            assertTrue(true);
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
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        return params;     
    }   
}
