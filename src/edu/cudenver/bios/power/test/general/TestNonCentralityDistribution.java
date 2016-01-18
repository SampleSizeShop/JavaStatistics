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
package edu.cudenver.bios.power.test.general;

import java.text.DecimalFormat;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.RandomColumnMetaData;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.NonCentralityDistribution;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import jsc.distributions.Normal;
import junit.framework.TestCase;

/**
 * Non-automated test of non-centrality distribution.  Manually compared against
 * Glueck & Muller 2003 implementation in SAS/IML
 * @author Sarah Kreidler
 *
 */
public class TestNonCentralityDistribution extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 1.0;
    private static final double[] ALPHA_LIST = {0.05};
    private static final double[] BETA_SCALE_LIST = {1};
    private static final double[] SIGMA_SCALE_LIST = {1};
    private static final int[] SAMPLE_SIZE_LIST = {5};
    private Normal normalDist = new Normal();
    private DecimalFormat Number = new DecimalFormat("#0.0000000000");

    public void testApproximateNonCentralCDF()
    {
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        double betaScale = params.getBetaScaleList().get(0);
        double sigmaScale = params.getSigmaScaleList().get(0);
        int perGroupN = params.getSampleSizeList().get(0);
        Test test = params.getTestList().get(0);

        //        Test test, RealMatrix F, RealMatrix FtFinverse, int N,
        //        FixedRandomMatrix CFixedRand, RealMatrix U,
        //        RealMatrix thetaNull, RealMatrix beta,
        //        RealMatrix sigmaError, RealMatrix sigmaG, boolean exact

        try {
            NonCentralityDistribution ncd =
                    new NonCentralityDistribution(test, params.getDesignEssence(),
                            null,
                            perGroupN,
                            params.getBetweenSubjectContrast(),
                            params.getWithinSubjectContrast(),
                            params.getTheta(),
                            params.getBeta().scalarMultiply(betaScale, true),
                            params.getSigmaError().scalarMultiply(sigmaScale),
                            params.getSigmaGaussianRandom(),
                            false);

            for(double crit = 1; crit < 15; crit++)
            {
                double prob = ncd.cdf(crit);
                System.out.println("Critical Value: " + crit + " prob: " + Number.format(prob));
            }
        } catch (PowerException e) {
            System.err.println("CDF failed: [" + e.getErrorCode() +
                    "]" + e.getMessage());
        }
    }

    public void testApproximateNonCentralInverseCDF()
    {
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        double betaScale = params.getBetaScaleList().get(0);
        double sigmaScale = params.getSigmaScaleList().get(0);
        int perGroupN = params.getSampleSizeList().get(0);
        Test test = params.getTestList().get(0);
        try {
            NonCentralityDistribution ncd =
                    new NonCentralityDistribution(test, params.getDesignEssence(),
                            null,
                            perGroupN,
                            params.getBetweenSubjectContrast(),
                            params.getWithinSubjectContrast(),
                            params.getTheta(),
                            params.getBeta().scalarMultiply(betaScale, true),
                            params.getSigmaError().scalarMultiply(sigmaScale),
                            params.getSigmaGaussianRandom(),
                            false);

            for(double w = 0.10; w < 1.05; w += 0.1)
            {
                double nonCentralityParam = ncd.inverseCDF(w);
                System.out.println("Quantile: " + Number.format(w) +
                        " inverseCDF: " + Number.format(nonCentralityParam));
            }
        } catch (PowerException e) {
            System.err.println("Simulation failed: [" + e.getErrorCode() +
                    "]" + e.getMessage());
        }
    }

    /**
     * Builds matrices for a multivariate GLM with a baseline covariate
     */
    private GLMMPowerParameters buildValidMultivariateRandomInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();

        // add power methods
        //for(PowerMethod method: PowerMethod.values()) params.addPowerMethod(method);
        params.addPowerMethod(PowerMethod.CONDITIONAL_POWER);
        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
        params.addQuantile(0.25);
        params.addQuantile(0.5);
        params.addQuantile(0.75);

        // add tests - only HL andUNIREP value for random case
        params.addTest(Test.HOTELLING_LAWLEY_TRACE);
        //params.addTest(GLMMPowerParameters.Test.UNIREP);

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int P = 3;
        int Q = 3;
        // create design matrix
        params.setDesignEssence(org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(Q));
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);

        // build sigma G matrix
        double[][] sigmaG = {{1}};
        params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

        // build sigma Y matrix
        double rho = 0.4;
        double [][] sigmaY = {{1,0},{0,1}};
        params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

        // build sigma YG
        double rhoYG = 0.8;
        double [][] sigmaYG = {{0.9},{0}};
        params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);

        // build beta matrix
        double [][] beta = {{1,0},{0,0},{0,0}};
        double [][] betaRandom = {{0.9,0}};
        params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);

        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0},{1,0,-1}};
        double[][] betweenRandom = {{1},{1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));

        // build within subject contrast
        double [][] within = {{1,0},{0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        // set the sigma error matrix to [sigmaY - sigmaYG * sigmaG-1 * sigmaGY]
        RealMatrix sigmaGY = params.getSigmaOutcomeGaussianRandom().transpose();
        RealMatrix sigmaGInverse = new LUDecomposition(params.getSigmaGaussianRandom()).getSolver().getInverse();
        params.setSigmaError(params.getSigmaOutcome().subtract(params.getSigmaOutcomeGaussianRandom().multiply(sigmaGInverse.multiply(sigmaGY))));

        return params;
    }
}
