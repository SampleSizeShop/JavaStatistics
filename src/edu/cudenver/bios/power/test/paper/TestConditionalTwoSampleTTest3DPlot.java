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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import junit.framework.TestCase;

/**
 * Unit test for fixed univariate design with comparison against
 * simulation and SAS output.
 * 
 * based on the example 3 from POWERLIB:
*   Johnson J.L., Muller K.E., Slaughter J.C., Gurka M.J., Gribbin M.J. and Simpson S.L. 
*   (2009) POWERLIB: SAS/IML software for computing power in multivariate linear models, 
*   Journal of Statistical Software, 30(5), 1-27.
*   
 * @author Sarah Kreidler
 *
 */
public class TestConditionalTwoSampleTTest3DPlot extends TestCase
{
	private static final String DATA_FILE =  "data" + File.separator + "TestConditionalTwoSampleTTest3DPlot.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + 
	        File.separator + "TestConditionalTwoSampleTTest3DPlot.tex";
	private static final String TITLE = "GLMM(F) Example 3. Power for a two sample " +
			"t-test for various sample sizes and mean differences";
	private static final String IMAGE_FILE = "TestConditionalTwoSampleTTest3DPlot.png";
	private static final String AUTHOR = "Sarah Kreidler";
	private static final String STUDY_DESIGN_DESCRIPTION  = 
	        "The study design for Example 3 is a balanced, two sample design with" +
            "a single response variable. We calculate power for a two-sample t-test comparing " +
            "the mean responses between the two independent groups.  " +
            "The example demonstrates changes in power with different " +
            "sample sizes and mean differences.";
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
     * Compare 2 sample t-test results between JavaStatistics, 
     * POWERLIB, and simulation
     */
    public void testTwoSampleTTEst()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests
        params.addTest(Test.UNIREP);
        
        // add alpha values
        params.addAlpha(0.01);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        for(double betaScale = 0; betaScale <= 0.75; betaScale += 0.05) params.addBetaScale(betaScale);
        
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
        for(int size = 3; size <= 18; size += 3) params.addSampleSize(size);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        checker.checkPower(params);
        // output the results
        try {
            ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
            reportBuilder.createValidationReportAsStdout(checker, TITLE, false);
            reportBuilder.createValidationReportAsLaTex(
                    OUTPUT_FILE, TITLE, AUTHOR, STUDY_DESIGN_DESCRIPTION, 
                    params, checker, IMAGE_FILE);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

		assertTrue(checker.isSASDeviationBelowTolerance());
		checker.reset();
    }
}
