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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;

/**
 * Test case for exact quantile power for the UNIREP.  Values should match
 * exact median power values from Table II in Glueck & Muller 2003.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestUnirepExactQuantile extends TestCase 
{
	private static final double[] ALPHA_LIST = {0.05};    
	private static final double[] SIGMA_SCALE_LIST = {1};	


	private static final String DATA_FILE =  "data" + File.separator + "TestUnirepExactQuantile.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + 
	File.separator + "UnirepExactQuantileOutput.tex";
	private static final String TITLE = "GLMM(F, g) Example 6. Median power for the uncorrected " +
			"univariate approach to repeated measures, Box, Geisser-Greenhouse," +
			" and Huynh-Feldt tests, using Davies Algorithm";
    private static final String AUTHOR = "Sarah Kreidler";
    private static final String STUDY_DESIGN_DESCRIPTION  = 
            "The study design in Example 6 is a three sample design with " +
                    "a baseline covariate and four repeated measurements.  We calculate " +
                    "the median power for a test of no difference between groups at each " +
                    "time point.  We calculate median power for the " +
                    "uncorrected univariate approach to repeated measures, " +
                    "Box, Geisser-Greenhouse, and Huynh-Feldt tests." +
                    "The exact distribution of the test statistic under the alternative hypothesis is obtained " +
                    "using Davies' algorithm described in \n\n" +
                    "\\hangindent2em\n\\hangafter=1\n Davies, R. B. (1980). " +
                    "Algorithm AS 155: The Distribution of a Linear Combination " +
                    "of Chi-Square Random Variables. \\emph{Applied Statistics}, " +
                    "\\emph{29}(3), 323-333.\n\n" +
                    "Median power is calculated for the following combinations " +
                    "of mean differences and per group sample sizes.\n\n" +
                    "\\begin{enumerate}" +
                    "\\item Per group sample size of 5, with beta scale values " +
                    "0.4997025, 0.8075886, and 1.097641" +
                    "\\item Per group sample size of 25, with beta scale values " +
                    "0.1651525, 0.2623301, and 0.3508015" +
                    "\\item Per group sample size of 50, with beta scale values " +
                    "0.1141548,  0.1812892, and  0.2423835\n" +
                    "\\end{enumerate}\n\n" +
                    "The example is based on Table II from\n\n" +
                    "\\hangindent2em\n\\hangafter=1\n Glueck, D. H., \\& Muller, K. E. (2003). " +
                    "Adjusting power for a baseline covariate in linear models. \\emph{Statistics " +
                    "in Medicine}, \\emph{22}(16), 2535-2551.\n\n";
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

		Test[] list = {Test.UNIREP,Test.UNIREP_BOX,
				Test.UNIREP_GEISSER_GREENHOUSE,Test.UNIREP_HUYNH_FELDT};
		for(Test test : list)
		{
			params5.clearTestList();
			params5.addTest(test);
			params25.clearTestList();
			params25.addTest(test);
			params50.clearTestList();
			params50.addTest(test);
			checker.checkPower(params5);
			checker.checkPower(params25);
			checker.checkPower(params50);
		}
        // output the results
        try {
            // clear the beta scale list and per group N since this is described in the
            // study design section and may be confusing if we list all the beta scales
            // twice.
            params50.clearBetaScaleList();
            params50.clearSampleSizeList();
            ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
            reportBuilder.createValidationReportAsStdout(checker, TITLE, true);
            reportBuilder.createValidationReportAsLaTex(
                    OUTPUT_FILE, TITLE, AUTHOR, STUDY_DESIGN_DESCRIPTION, 
                    params50, checker);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
		assertTrue(checker.isSASDeviationBelowTolerance());
		checker.reset();	
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

		// add alpha values
		for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

		int P = 3;
		int Q = 3;
		// create design matrix
		params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
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
