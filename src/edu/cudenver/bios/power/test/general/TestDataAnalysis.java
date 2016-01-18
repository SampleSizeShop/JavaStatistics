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
import org.apache.commons.math3.linear.SingularValueDecomposition;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateEpsilonApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.ModelFit;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.glmm.GLMMTestHotellingLawley;
import edu.cudenver.bios.power.glmm.GLMMTestPillaiBartlett;
import edu.cudenver.bios.power.glmm.GLMMTestUnirepBox;
import edu.cudenver.bios.power.glmm.GLMMTestUnirepGeisserGreenhouse;
import edu.cudenver.bios.power.glmm.GLMMTestUnirepHuynhFeldt;
import edu.cudenver.bios.power.glmm.GLMMTestUnivariateRepeatedMeasures;
import edu.cudenver.bios.power.glmm.GLMMTestWilksLambda;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import jsc.distributions.FishersF;
import jsc.distributions.Normal;
import junit.framework.TestCase;

/**
 * Test the GLMMTest objects observed F, degrees of freedom when
 * doing data analysis - used for simulation.  This isn't really an
 * automated unit test - used to compare manually against LINMOD
 * library (SAS/IML) written by Keith Muller.
 *
 * @author Sarah Kreidler
 *
 */
public class TestDataAnalysis extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};
    private static final double[] BETA_SCALE_LIST = {2};
    private static final double[] SIGMA_SCALE_LIST = {2};
    private static final int[] SAMPLE_SIZE_LIST = {5};
    private Normal normalDist = new Normal();
    private DecimalFormat number = new DecimalFormat("#0.000");

    // design matrix for multivariate tests
    private RealMatrix essenceMultivariate =
        org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(4);
    private RealMatrix XMultivariate = MatrixUtils.getKroneckerProduct(essenceMultivariate,
            MatrixUtils.getRealMatrixWithFilledValue(5, 1, 1));
    private int rankXMultivariate = new SingularValueDecomposition(essenceMultivariate).getRank();
    private RealMatrix XtXInverseMultivariate =
        new LUDecomposition(XMultivariate.transpose().multiply(XMultivariate)).getSolver().getInverse();

    // design matrix for univariate tests
    private RealMatrix essenceUnivariate =
        org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(2);
    private RealMatrix XUnivariate = MatrixUtils.getKroneckerProduct(essenceUnivariate,
            MatrixUtils.getRealMatrixWithFilledValue(10, 1, 1));
    private int rankXUnivariate = new SingularValueDecomposition(essenceUnivariate).getRank();
    private RealMatrix XtXInverseUnivariate =
        new LUDecomposition(XUnivariate.transpose().multiply(XUnivariate)).getSolver().getInverse();

    // contrasts - univariate
    double [][] betweenUnivariate = {{1,-1}};
    RealMatrix CUnivariate = new Array2DRowRealMatrix(betweenUnivariate);
    RealMatrix UUnivariate = MatrixUtils.getRealMatrixWithFilledValue(1, 1, 1);
    RealMatrix ThetaNullUnivariate = MatrixUtils.getRealMatrixWithFilledValue(1, 1, 0);
    // contrasts - multivariate
    double [][] betweenMultivariate = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
    RealMatrix CMultivariate = new Array2DRowRealMatrix(betweenMultivariate);
    double [][] withinMultivariate = {{1,1},{-1,0},{0,-1}};
    RealMatrix UMultivariate = new Array2DRowRealMatrix(withinMultivariate);
    RealMatrix ThetaNullMultivariate = MatrixUtils.getRealMatrixWithFilledValue(3, 2, 0);

    // univariate Y data
    double[] YUnivariateData =
    {
            2.67,
            0.97,
            0.55,
            1.57,
            -0.12,
            -1.02,
            0.94,
            -1.25,
            0.94,
            -0.08,
            2.09,
            2.7800000000000002,
            1.7,
            1.4100000000000001,
            1.47,
            2.2,
            4.890000000000001,
            0.96,
            4.05,
            0.9199999999999999
    };
    private RealMatrix YUnivariate = new Array2DRowRealMatrix(YUnivariateData);

    // multivariate Y data
    double[][] YMultivariateData =
    {
            {0.43999999999999995, -0.88, 0.58},
            {4.029999999999999, -1.88, -0.16},
            {2.19, -0.07, 0.67},
            {2.83, -0.22, 0.63},
            {0.6000000000000001, -1.47, 1.71},
            {0.07, -0.53, -0.09},
            {-0.2, -0.87, 0.06},
            {0.27, -1.12, -1.42},
            {-0.01, 0.08, -0.38},
            {1.66, -0.53, -0.75},
            {0.25, 1.73, 0.08},
            {-1.05, 0.61, -0.15},
            {2.01, -0.54, -1.34},
            {-1.06, 0.16, 1.88},
            {-1.43, 2.18, -0.39},
            {-0.75, -0.19, -1.46},
            {-1.08, -0.91, 0.8},
            {-1.16, -1.68, 0.61},
            {1.67, 1.0, 0.55},
            {-0.47, -0.07, 0.14}
    };
    private RealMatrix YMultivariate = new Array2DRowRealMatrix(YMultivariateData);

    public void testUnirep()
    {
        GLMMTest unirepTest = new GLMMTestUnivariateRepeatedMeasures(
                FApproximation.NONE,
                UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                XUnivariate, XtXInverseUnivariate, rankXUnivariate, YUnivariate,
        CUnivariate, UUnivariate, ThetaNullUnivariate);

        ModelFit fit = unirepTest.getModelFit();
        checkFit(Test.UNIREP, fit);
    }

    public void testBox()
    {
        GLMMTest unirepBoxTest = new GLMMTestUnirepBox(
                FApproximation.NONE,
                UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                XUnivariate, XtXInverseUnivariate, rankXUnivariate, YUnivariate,
        CUnivariate, UUnivariate, ThetaNullUnivariate);

        ModelFit fit = unirepBoxTest.getModelFit();
        checkFit(Test.UNIREP_BOX, fit);

    }

    public void testGeisserGreenhouse()
    {
        GLMMTest unirepGGTest = new GLMMTestUnirepGeisserGreenhouse(
                FApproximation.NONE,
                UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                XUnivariate, XtXInverseUnivariate, rankXUnivariate, YUnivariate,
        CUnivariate, UUnivariate, ThetaNullUnivariate);

        ModelFit fit = unirepGGTest.getModelFit();
        checkFit(Test.UNIREP_GEISSER_GREENHOUSE, fit);

    }

    public void testHuynhFeldt()
    {
        GLMMTest unirepHFTest = new GLMMTestUnirepHuynhFeldt(
                FApproximation.NONE,
                UnivariateCdfApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                UnivariateEpsilonApproximation.MULLER_EDWARDS_TAYLOR_APPROX,
                XUnivariate, XtXInverseUnivariate, rankXUnivariate, YUnivariate,
        CUnivariate, UUnivariate, ThetaNullUnivariate);

        ModelFit fit = unirepHFTest.getModelFit();
        checkFit(Test.UNIREP_HUYNH_FELDT, fit);
    }

    public void testWilksLambda()
    {
        GLMMTest wilksTest = new GLMMTestWilksLambda(
                FApproximation.NONE,
                XMultivariate, XtXInverseMultivariate, rankXMultivariate, YMultivariate,
        CMultivariate, UMultivariate, ThetaNullMultivariate);

        ModelFit fit = wilksTest.getModelFit();
        checkFit(Test.WILKS_LAMBDA, fit);
    }

    public void testPillaiBartlett()
    {
        GLMMTest pillaiTest = new GLMMTestPillaiBartlett(
                FApproximation.NONE,
                XMultivariate, XtXInverseMultivariate, rankXMultivariate, YMultivariate,
        CMultivariate, UMultivariate, ThetaNullMultivariate);

        ModelFit fit = pillaiTest.getModelFit();
        checkFit(Test.PILLAI_BARTLETT_TRACE, fit);
    }

    public void testHotellingLawley()
    {
        GLMMTest hltTest = new GLMMTestHotellingLawley(
                FApproximation.NONE,
                XMultivariate, XtXInverseMultivariate, rankXMultivariate, YMultivariate,
        CMultivariate, UMultivariate, ThetaNullMultivariate);

        ModelFit fit = hltTest.getModelFit();
        checkFit(Test.HOTELLING_LAWLEY_TRACE, fit);
    }

    private void checkFit(Test test, ModelFit fit)
    {
        System.out.println("Test: " + test +
                " Ndf: " + number.format(fit.numeratorDF) +
                " Ddf: " + number.format(fit.denominatorDF) +
                " F-crit: " + number.format(fit.Fvalue) +
                " p-value: " + number.format(fit.Pvalue));
        RealMatrix beta = fit.beta;
        System.out.println("Parameter estimates:");
        for(int r = 0; r < beta.getRowDimension(); r++)
        {
            for(int c = 0; c < beta.getColumnDimension(); c++)
            {
                System.out.println("Beta["+ r+ "," + c+ "]=" + beta.getEntry(r, c));
            }
        }
    }
}
