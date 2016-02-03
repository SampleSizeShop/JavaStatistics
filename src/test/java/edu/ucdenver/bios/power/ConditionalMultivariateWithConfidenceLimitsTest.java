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
import edu.cudenver.bios.matrix.OrthogonalPolynomials;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import edu.cudenver.bios.utils.Factor;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * Unit test for fixed multivariate design including confidence intervals
 * with comparison against simulation and SAS output.
 *
 *  based on the example 6 from POWERLIB:
*   Johnson J.L., Muller K.E., Slaughter J.C., Gurka M.J., Gribbin M.J. and Simpson S.L.
*   (2009) POWERLIB: SAS/IML software for computing power in multivariate linear models,
*   Journal of Statistical Software, 30(5), 1-27.
*
*   !WARNING! - you will need to increase the main_memory size in your
*   LaTeX install to typeset this report.  To do so:
*   1. Open a terminal window
*   2. At the prompt enter
*      initexmf --edit-config-file=latex
*   3. In the editor window, type
*   main_memory=10000000
*   Save the file then quit the editor
*   4. Rebuild the format by running
*   initexmf --dump=latex
*   5. Repeat steps 2-4 with config files 'pdflatex' and 'xelatex'
*
 * @author Sarah Kreidler
 *
 */
public class ConditionalMultivariateWithConfidenceLimitsTest {

    private static final String DATA_FILE =  "TestConditionalMultivariateWithConfidenceLimits.xml";
    private static final String TITLE = "GLMM(F) Example 6. Power and confidence " +
            "limits for the univariate approach to repeated measures in a multivariate model";
    private static final double TOLERANCE = 0.000001;

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

    @Before
    public void setUp() {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, false);
    }

    /**
     * Test valid inputs for a univariate linear model with only fixed predictors
     */
    @Test
    public void testMultivariateWithConfidenceLimits()
    {
        GLMMTestFactory.Test[] testList = {GLMMTestFactory.Test.UNIREP_GEISSER_GREENHOUSE};
        for(GLMMTestFactory.Test test: testList)
        {
            GLMMPowerParameters params = buildInputs(test);

            for(double delta = 0.0008; delta < 0.2001; delta += 0.0008)
            {
                // increase the gender difference by 2 * delta
                RealMatrix betaMatrix = params.getBeta().getFixedMatrix();
                for(int row = 0; row < 5; row++) betaMatrix.setEntry(row, 2, beta[row][2] + delta);
                for(int row = 5; row < 10; row++) betaMatrix.setEntry(row, 2, beta[row][2] - delta);

                checker.checkPower(params);
            }
        }

        // we add all of the tests back into the params to ensure the report lists them all
        // we do this to ensure the order of results in the SAS output matches the order
        // in the JavaStatistics calls.
        GLMMPowerParameters params = buildInputs(GLMMTestFactory.Test.UNIREP);
        params.clearTestList();
        for(GLMMTestFactory.Test test: testList) {
            params.addTest(test);
        }

        // output the results
        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, false);
        assertTrue("SAS deviation " + checker.getMaxSasDeviation() + " is not below tolerance " + TOLERANCE,
                    checker.isSASDeviationBelowTolerance(TOLERANCE));
    }

    private GLMMPowerParameters buildInputs(GLMMTestFactory.Test test)
    {
        // build the inputs
        GLMMPowerParameters params = new GLMMPowerParameters();

        // add tests
        params.addTest(test);

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
        params.setDesignEssence(org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(10));
        // add sample size multipliers
        for(int sampleSize = 2; sampleSize <= 10; sampleSize++) params.addSampleSize(sampleSize);

        // build beta matrix
        params.setBeta(new FixedRandomMatrix(beta, null, false));

        // build between subject contrast
        double [][] between = {{1,1,1,1,1,-1,-1,-1,-1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        double[] regions = {1,2,3,4};
        String name = "region";
        ArrayList<Factor> factorList = new ArrayList<Factor>();
        Factor regionFactor = new Factor(name, regions);
        factorList.add(regionFactor);
        params.setWithinSubjectContrast(
                OrthogonalPolynomials.withinSubjectContrast(factorList).getMainEffectContrast(regionFactor).getContrastMatrix());

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
