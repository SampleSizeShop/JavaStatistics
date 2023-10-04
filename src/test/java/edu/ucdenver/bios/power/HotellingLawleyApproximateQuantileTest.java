/*
 * Java Statistics.  A java library providing power/sample size estimation for
 * the general linear model.
 *
 * Copyright (C) 2016 Regents of the University of Colorado.
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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;

import javax.xml.parsers.ParserConfigurationException;

import static org.junit.Assert.assertTrue;

/**
 * Test case for approximate quantile power for the HLT.  Values should match
 * approximate median power values from Table II in Glueck & Muller 2003.
 *
 * @author Sarah Kreidler
 *
 */
public class HotellingLawleyApproximateQuantileTest {
    private static final String DATA_FILE =  "TestHotellingLawleyApproximateQuantile.xml";
    private static final String OUTPUT_FILE = "text" + File.separator + "results" +
    File.separator + "HotellingLawleyApproximateQuantileOutput.tex";
    private static final String TITLE = "GLMM(F, g) Example 1. Median power for the " +
            "Hotelling-Lawley Trace, using the Satterthwaite approximation";
    private static final double[] ALPHA_LIST = {0.05};
    private static final double[] SIGMA_SCALE_LIST = {1};
    private static final double TOLERANCE = 5e-5;

    private PowerChecker checker;

    @Before
    public void setUp() throws ParserConfigurationException, SAXException {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, true);
    }

    /**
     * Compare the calculated HLT approximate quantile powers against simulation
     */
    @Test
    public void testPower() {
        // build the inputs
        double[] beta5 = {
                0.4997025,
                0.8075886,
                1.097641};
        GLMMPowerParameters params5 = validMultivariateRandomInputs(beta5, 5);
        double[] beta25 = {
                0.1651525,
                0.2623301,
                0.3508015
        };
        GLMMPowerParameters params25 = validMultivariateRandomInputs(beta25, 25);
        double[] beta50 = {
                0.1141548,
                0.1812892,
                0.2423835
        };
        GLMMPowerParameters params50 = validMultivariateRandomInputs(beta50, 50);

        checker.checkPower(params5);
        checker.checkPower(params25);
        checker.checkPower(params50);

        // output the results
        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, true);

        assertTrue("results outside tolerance: " + TOLERANCE, checker.isSASDeviationBelowTolerance(TOLERANCE));
        checker.reset();
    }

    /**
     * Builds matrices for a multivariate GLM with a baseline covariate.
     *
     * <p>
     * Note: this matrix set matches the values produced in Table II from Glueck&Muller.
     *
     * @return The power parameters reflecting the matrices.
     */
    private GLMMPowerParameters validMultivariateRandomInputs() {
        return validMultivariateRandomInputs(null, 0);
    }

    /**
     * Builds matrices for a multivariate GLM with a baseline covariate, incorporating
     * a beta scale list and a sample size.
     *
     * <p>
     * Note: this matrix set matches the values produced in Table II from Glueck&Muller.
     *
     * @param betaScaleList The beta scale list.
     * @param repn          The sample size.
     *
     * @return The power parameters reflecting the matrices, the beta scale list, and the
     *         sample size.
     */
    private GLMMPowerParameters validMultivariateRandomInputs(double[] betaScaleList, int repn) {
        GLMMPowerParameters params = new GLMMPowerParameters();

        // add quantile power method
        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
        params.addQuantile(0.5);

        // add HLT as the statistical test
        params.addTest(GLMMTestFactory.Test.HOTELLING_LAWLEY_TRACE);

        // add alpha values
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        // create design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(3));

        // add sample size multipliers
        //  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        if (repn != 0) {
            params.addSampleSize(repn);
        }

        // build sigma G matrix
        double[][] sigmaG = {{1}};
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
        if (betaScaleList != null) {
            for(double betaScale: betaScaleList) params.addBetaScale(betaScale);
        }

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
