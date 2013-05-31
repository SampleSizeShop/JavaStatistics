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

import junit.framework.TestCase;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;

/**
 * Unit test for fixed univariate design with confidence intervals (powerlib example 4). 
 * Compared against simulation and SAS output.
 * 
 *  based on the example 4 from POWERLIB:
*   Johnson J.L., Muller K.E., Slaughter J.C., Gurka M.J., Gribbin M.J. and Simpson S.L. 
*   (2009) POWERLIB: SAS/IML software for computing power in multivariate linear models, 
*   Journal of Statistical Software, 30(5), 1-27.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestConditionalUnivariateWithConfidenceLimits extends TestCase
{
	private static final String DATA_FILE = "data" + File.separator + 
	        "TestConditionalUnivariateWithConfidenceLimits.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + 
	        File.separator + "TestConditionalUnivariateWithConfidenceLimits.tex";
	private static final String TITLE = "GLMM(F) Example 4. Power and confidence limits for a univariate model";
	   private static final String AUTHOR = "Sarah Kreidler";
	    private static final String STUDY_DESIGN_DESCRIPTION  = 
	            "The study design for Example 4 is a balanced two group design. " +
	            "We calculate power for a two-sample t-test comparing the mean response " +
	            "between the groups.  We calculate confidence limits for the power values. " +
	            "The example is based on Figure 1 " +
	            "from\n\n\\hangindent2em\n\\hangafter=1\n Taylor, D. J., \\& Muller, K. E. (1995). " +
	            "Computing Confidence Bounds for Power and Sample Size of the " +
	            "General Linear Univariate Model. \\emph{The American Statistician}, " +
	            "\\emph{49}(1), 43-47.";
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
			e.printStackTrace();
			fail();
		}
	}

    /**
     * Compare 2 sample t-test results between JavaStatistics, 
     * POWERLIB, and simulation
     */
    public void testUnviariateWithConfidenceLimits()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests
        params.addTest(Test.UNIREP);
        
        // add alpha values
        params.addAlpha(0.01);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values from 0 to 0.75
        for(double betaScale = 0; betaScale < 0.76; betaScale += 0.01) params.addBetaScale(betaScale);

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
        params.setConfidenceIntervalType(ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED);
        params.setSampleSizeForEstimates(24);
        params.setDesignMatrixRankForEstimates(2);
        
        // run the test
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
        // output the results
        try {
            ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
            reportBuilder.createValidationReportAsStdout(checker, TITLE, false);
            reportBuilder.createValidationReportAsLaTex(
                    OUTPUT_FILE, TITLE, AUTHOR, STUDY_DESIGN_DESCRIPTION, 
                    params, checker);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

		assertTrue(checker.isSASDeviationBelowTolerance());
		checker.reset();
    }

}
