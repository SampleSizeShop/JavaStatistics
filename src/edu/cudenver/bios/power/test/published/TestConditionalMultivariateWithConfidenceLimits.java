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

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.utils.OrthogonalPolynomials;
import junit.framework.TestCase;

/**
 * Unit test for fixed multivariate design including confidence intervals
 * with comparison against simulation and SAS output.
 * 
 *  based on the example 6 from POWERLIB:
*   Johnson J.L., Muller K.E., Slaughter J.C., Gurka M.J., Gribbin M.J. and Simpson S.L. 
*   (2009) POWERLIB: SAS/IML software for computing power in multivariate linear models, 
*   Journal of Statistical Software, 30(5), 1-27.
*   
 * @author Sarah Kreidler
 *
 */
public class TestConditionalMultivariateWithConfidenceLimits extends TestCase
{
	private static final String DATA_FILE =  "sas" + File.separator + "data" + File.separator + "TestConditionalMultivariateWithConfidenceLimits.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestConditionalMultivariateWithConfidenceLimits.html";
	private static final String TITLE = "Power results for multivariate design, UNIREP-GG with confidence limits";
	private PowerChecker checker;
	
    // set beta matrix
	private double[][] beta = 
	{
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2},
			{2.9, 3.2, 3.5, 3.2}
	};
	
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
     * Test valid inputs for a univariate linear model with only fixed predictors
     */
    public void testMultivariateWithConfidenceLimits()
    {
    	GLMMPowerParameters params = buildInputs();
    	
    	for(double delta = 0; delta < 0.0501; delta += 0.0008)
    	{
    		// increase the gender difference by 2 * delta
    		RealMatrix betaMatrix = params.getBeta().getFixedMatrix();
    		for(int row = 0; row < 5; row++) betaMatrix.setEntry(row, 2, beta[row][2] + delta);
    		for(int row = 5; row < 10; row++) betaMatrix.setEntry(row, 2, beta[row][2] - delta);

            checker.checkPower(params);
    	}

        System.out.println(TITLE);
		checker.outputResults();
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
		assertTrue(checker.isSimulationDeviationBelowTolerance());
		checker.reset();
    }

    
    private GLMMPowerParameters buildInputs()
    {
        // build the inputs        
    	GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests
    	//params.addTest(Test.UNIREP_GEISSER_GREENHOUSE);
    	params.addTest(Test.UNIREP);
    	params.addTest(Test.UNIREP_BOX);
    	params.addTest(Test.UNIREP_GEISSER_GREENHOUSE);
    	params.addTest(Test.UNIREP_HUYNH_FELDT);
    	params.addTest(Test.WILKS_LAMBDA);
    	params.addTest(Test.PILLAI_BARTLETT_TRACE);
    	params.addTest(Test.HOTELLING_LAWLEY_TRACE);

    	//params.addTest(Test.HOTELLING_LAWLEY_TRACE);
        // add alpha values - bonferroni corrected for 6 comparisons
        params.addAlpha(0.05/6);
        
        // add beta scale values
        params.addBetaScale(1);
        
        // build theta null matrix
        double [][] theta0 = {{0,0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double [][] sigma = {{0.08380, 0.05020, 0.03560, 0.05330},
       		 {0.05020, 0.05370, 0.03250, 0.03330},                         
             {0.03560, 0.03250, 0.04410, 0.03860},                          
             {0.05330, 0.03330, 0.03860, 0.07220}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        params.addSigmaScale(1);
        
        // build design matrix
        params.setDesignEssence(org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(10));
        // add sample size multipliers
        for(int sampleSize = 2; sampleSize <= 10; sampleSize++) params.addSampleSize(sampleSize);
        
        // build beta matrix
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        
        // build between subject contrast
        double [][] between = {{1,1,1,1,1,-1,-1,-1,-1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));
        
        double[] regions = {1,2,3,4};
        params.setWithinSubjectContrast(OrthogonalPolynomials.withinSubjectContrast(regions));
        
        // parameters for confidence limits
        params.setConfidenceIntervalType(ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED);
        params.setSampleSizeForEstimates(21);
        params.setDesignMatrixRankForEstimates(1);
        // 2 sided CI
        params.setAlphaLowerConfidenceLimit(0.025);
        params.setAlphaUpperConfidenceLimit(0.025);
        
        return params;
    }
}
