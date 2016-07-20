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

import java.io.File;
import java.util.List;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import junit.framework.TestCase;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Test case for exact unconditional power for the HLT.  Values should match
 * exact unconditional power values from Table II in Glueck & Muller 2003.
 *
 * @author Sarah Kreidler
 *
 */
public class HotellingLawleyExactUnconditionalTest {

    private static final String DATA_FILE =  "TestHotellingLawleyExactUnconditional.xml";
    private static final String TITLE = "GLMM(F, g) Example 4. Unconditional power for the " +
            "Hotelling-Lawley Trace, using Davies algorithm";
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 1.0;
    private static final double[] ALPHA_LIST = {0.05};
    private static final double[] SIGMA_SCALE_LIST = {1};
    private static final double TOLERANCE = 0.001;

    private PowerChecker checker;

    @Before
    public void setUp() {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, false);
    }
    /**
     * Compare the calculated HLT exact unconditional powers against simulation
     */
    @Test
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

        checker.checkPower(params5);
        checker.checkPower(params25);
        checker.checkPower(params50);

        // clear the beta scale list and per group N since this is described in the
        // study design section and may be confusing if we list all the beta scales
        // twice.
        params50.clearBetaScaleList();
        params50.clearSampleSizeList();
        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, false);
        assertTrue("results outside tolerance: " + TOLERANCE, checker.isSASDeviationBelowTolerance(TOLERANCE));
    }

    /**
     * Builds matrices for a multivariate GLM with a baseline covariate
     * Note, this matrix set matches the values produced in Table II from Glueck&Muller
     */
    private GLMMPowerParameters buildValidMultivariateRandomInputs(double[] betaScaleList, int repn)
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        params.setNonCentralityCDFExact(true);

        // add unconditional power methods and median unconditional
        params.addPowerMethod(PowerMethod.UNCONDITIONAL_POWER);

        // add HLT as the statistical test
        params.addTest(GLMMTestFactory.Test.HOTELLING_LAWLEY_TRACE);

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        // create design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(3));
        // add sample size multipliers
      //  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        params.addSampleSize(repn);
        // build sigma G matrix
        double[][] sigmaG = {{VARIANCE}};
        params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

        // build sigma Y matrix
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
        double [][] between = {{-1,1,0}, {-1,0,1}};
        double[][] betweenRandom = {{0}, {0}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));

        // build within subject contrast
        double [][] within = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;
    }
}
