/*
 * Java Statistics.  A java library providing power/sample size estimation for
 * the general linear model.
 *
 * Copyright (C) 2017 Regents of the University of Colorado.
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
package edu.cudenver.bios.power;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.StatUtils;

import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtilities;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.RandomErrorMatrix;
import edu.cudenver.bios.power.DetectableDifferenceBound.DetectableDifferenceError;
import edu.cudenver.bios.power.SampleSizeBound.SampleSizeError;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTest.ModelFit;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.glmm.NonCentralityDistribution;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.utils.Logger;

import static edu.cudenver.bios.matrix.MatrixUtilities.forceSymmetric;

import java.math.BigInteger;

/**
 * Power calculator implementation for the general linear multivariate model
 *
 * @see PowerCalculator
 * @author Sarah Kreidler
 *
 */
public class GLMMPowerCalculator implements PowerCalculator
{
    private BigInteger numberOfEvaluations;

    private static final String NO_MEAN_DIFFERENCE_MESSAGE =
            "The null hypothesis is true: that is, all contrasts defined by the hypothesis have zero sums of squares. "
        +   "(This may arise, for example, in a test of mean difference if the means are equal.) "
        +   "Thus the highest possible power is \u03B1 (alpha, the Type I error rate), "
        +   "and no sample size can be large enough to achieve higher power."
        ;

    private static final String SIGMA_ERROR_NOT_POSITIVE_SEMIDEFINITE_MESSAGE =
            "Unfortunately, there is no solution for this combination of input parameters. "
        +   "The error covariance matrix (\u03a3<sub>E</sub>) does not describe a valid "
        +   "covariance structure (that is, it is not positive semidefinite). "
        +   "Reducing the expected covariate-to-response correlations "
        +   "will likely lead to a soluble combination."
        ;

    private static final String NEGATIVE_NU_EST_MESSAGE =
            "For confidence interval calculation, the total sample size must be greater than "
        +   "the rank of the design matrix. "
        +   "Please revisit the Options > Confidence Intervals page and make it so."
        ;

    private static final int MAX_ITERATIONS = 10000;

    private static final int MAX_EVALUATIONS = 0x200003; // 2^21 + 3: smallest needed for unit tests to pass

    private static final int STARTING_SAMPLE_SIZE = 0x400;     // 2^10
    private static final int MAX_SAMPLE_SIZE      = 0x2000000; // 2^25: largest for which performance is good

    private static final int STARTING_BETA_SCALE = 10;

    private static final int SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL = 1000;
    private static final double EPSILON = 1e-12;

    // The largest sample size for which we will agree to perform a power calculation
    // using the unconditional power method. (Beyond this, performance may suffer.)
    private static final int MAX_SAMPLE_SIZE_FOR_UNCONDITIONAL_POWER = 100;

    // seed for random column generation
    private int seed = 1234;
    // accuracy thresholds
    // minimum value still considered positive in Cholesky decomposition
    protected double positivityThreshold =
            CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
    // minimum difference still considered symmetric in Cholesky decomposition
    protected double symmetryThreshold =
            CholeskyDecomposition.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;

    private static final List<Double> NO_QUANTILES = Collections.singletonList(Double.NaN);

    private static final Logger LOGGER = Logger.getLogger(GLMMPowerCalculator.class);

    public class SimulatedPower
    {
        private double power;
        private RealMatrix averageBeta;
        private int iterations = 0;

        public SimulatedPower(int betaRow, int betaCol)
        {
            averageBeta = MatrixUtils.getRealMatrixWithFilledValue(betaRow, betaCol, 0);
        }

        public double getPower() { return power; }

        public void setPower(double power)    {this.power = power;}

        public RealMatrix getAverageBeta()
        {
            return averageBeta.scalarMultiply((double) 1/(double)this.iterations);
        }

        public void accumulateBeta(RealMatrix averageBeta)
        {
            this.averageBeta = this.averageBeta.add(averageBeta);
            iterations++;
        }
    }

    /**
     * container class for simulation info
     */
    private class SimulationFit
    {
        public RealMatrix Y;
        public RealMatrix Ypred;
        public RealMatrix sigma;
        public RealMatrix beta;
        public double numeratorDF;
        public double denominatorDF;
        public double Fvalue;
        public double Pvalue;

        public SimulationFit(double Pvalue, double Fvalue,
                double numeratorDF, double denominatorDF,
                RealMatrix Y, RealMatrix Ypred, RealMatrix sigma, RealMatrix beta)
        {
            this.Y = Y;
            this.sigma = sigma;
            this.beta = beta;
            this.Pvalue = Pvalue;
            this.Fvalue = Fvalue;
            this.numeratorDF = numeratorDF;
            this.denominatorDF = denominatorDF;
        }
    }

    /**
     * Function used with Apache's bisection solver to determine the
     * per-group sample size which most closely achieves the desired power
     */
    private class SampleSizeFunction implements UnivariateFunction
    {
        private GLMMTest glmmTest;
        private NonCentralityDistribution nonCentralityDist;
        private PowerMethod method;
        private double alpha;
        private double quantile;
        private double targetPower;

        public SampleSizeFunction(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
                PowerMethod method, double targetPower, double alpha, double quantile)
        {
            this.glmmTest = glmmTest;
            this.nonCentralityDist = nonCentralityDist;
            this.method = method;
            this.targetPower = targetPower;
            this.alpha = alpha;
            this.quantile = quantile;
        }

        public double value(double n)
        {
            try
            {
                glmmTest.setPerGroupSampleSize((int) n);
                if (nonCentralityDist != null) nonCentralityDist.setPerGroupSampleSize((int) n);
                double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
                double diff = targetPower - calculatedPower;
                return diff;
            }
            catch (Exception e)
            {
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a negative number
                return Double.NaN;
            }
        }
    }

    /**
     * Function used with Apache's bisection solver to determine the
     * per-group sample size which most closely achieves the desired power
     */
    private class DetectableDifferenceFunction implements UnivariateFunction
    {
        private GLMMTest glmmTest;
        private NonCentralityDistribution nonCentralityDist;
        private PowerMethod method;
        private double alpha;
        private double quantile;
        private double targetPower;
        private FixedRandomMatrix beta;

        public DetectableDifferenceFunction(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
                PowerMethod method, double targetPower, double alpha, double quantile, FixedRandomMatrix beta)
        {
            this.glmmTest = glmmTest;
            this.nonCentralityDist = nonCentralityDist;
            this.method = method;
            this.alpha = alpha;
            this.quantile = quantile;
            this.targetPower = targetPower;
            this.beta = beta;
        }

        public double value(double betaScale)
        {
            try
            {
                RealMatrix scaledBeta = beta.scalarMultiply(betaScale, true);
                glmmTest.setBeta(scaledBeta);
                if (nonCentralityDist != null) nonCentralityDist.setBeta(scaledBeta);
                double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
                double diff = targetPower - calculatedPower;
                return diff;
            }
            catch (Exception e)
            {
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a negative number
                return Double.NaN;
            }
        }
    }

    /**
     * Class passed into Apache's SimpsonIntegrator function to compute
     * unconditional power
     */
    private class UnconditionalPowerIntegrand implements UnivariateFunction
    {
        protected NonCentralityDistribution nonCentralityDist;
        protected double Fcrit;
        protected double ndf;
        protected double ddf;

        public UnconditionalPowerIntegrand(NonCentralityDistribution nonCentralityDist,
                double Fcrit, double ndf, double ddf)
        {
            this.nonCentralityDist = nonCentralityDist;
            this.Fcrit = Fcrit;
            this.ndf = ndf;
            this.ddf = ddf;
        }

        public double value(double t)
        {
            numberOfEvaluations = numberOfEvaluations.add(BigInteger.ONE);

            NonCentralFDistribution FdistTerm1 = new NonCentralFDistribution(ndf, ddf, t);
            NonCentralFDistribution FdistTerm2 = new NonCentralFDistribution(ndf+2, ddf, t);

            try {
                return nonCentralityDist.cdf(t)*(FdistTerm1.cdf(Fcrit) - FdistTerm2.cdf((Fcrit*ndf)/(ndf+2)));
            } catch (PowerException pe) {
                return Double.NaN;
            }
        }
    }

    /********* public methods for the power API ************/

    /**
     * Calculate a list of power values using the methodology described in
     * Muller & barton 1984, Muller LaVange Ramey & Ramey 1992,
     * Muller, Edwards, Simpson & Taylor 2007, and Glueck & Muller 2003.
     * Please visit samplesizeshop.org for a full list of references.
     * If the parameters contain lists of possible scale factors, statistical
     * tests, etc., then a power will be returned for each combination
     * of these factors.
     *
     * @see GLMMPowerParameters
     * @see GLMMPower
     * @param powerParams inputs to the power calculation
     * @return list of calculated power values.
     */
    @Override
    public List<Power> getPower(PowerParameters powerParams)
            throws PowerException {
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);

        // precalculate any computationally expensive matrices/constants,
        // update the parameters as needed - used for random covariates
        initialize(params);

        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();

        // calculate the power for all variations of the study design
        for (Test test: params.getTestList()) {
            for (PowerMethod method: params.getPowerMethodList()) {
                List<Double> quantileList = method == PowerMethod.QUANTILE_POWER
                                                ? params.getQuantileList()
                                                : NO_QUANTILES;

                for (Double alpha: params.getAlphaList()) {
                    for (Double sigmaScale: params.getSigmaScaleList()) {
                        for (Double betaScale: params.getBetaScaleList()) {
                            for (Double quantile: quantileList) {
                                for (Integer sampleSize: params.getSampleSizeList()) {
                                    /*
                                     * add the power result to the list
                                     * if a failure occurs, an error code and message are
                                     * included with this object
                                     */
                                    results.add(getPowerValue(params, test, method, alpha, sigmaScale, betaScale,
                                            sampleSize, quantile));
                                    if (Thread.currentThread().isInterrupted()) {
                                        return results;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return results;
    }

    /**
     * Find the best possible sample size to achieve a specified power or
     * list of powers.  Sample size is determined with a bisection search
     *
     * @see GLMMPowerParameters
     * @see GLMMPower
     * @param sampleSizeParams inputs to the sample size calculation
     * @return list of calculated power values.
     */
    @Override
    public List<Power> getSampleSize(PowerParameters sampleSizeParams)
            throws PowerException {
        GLMMPowerParameters params = (GLMMPowerParameters) sampleSizeParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);

        // precalculate any computationally expensive matrices/constants,
        // update the parameters as needed - used for random covariates
        initialize(params);

        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();

        // calculate the power for all variations of the study design
        for (Test test: params.getTestList()) {
            for (PowerMethod method: params.getPowerMethodList()) {
                List<Double> quantileList = method == PowerMethod.QUANTILE_POWER
                                                ? params.getQuantileList()
                                                : NO_QUANTILES;

                for (Double alpha: params.getAlphaList()) {
                    for (Double sigmaScale: params.getSigmaScaleList()) {
                        for (Double betaScale: params.getBetaScaleList()) {
                            for (Double quantile: quantileList) {
                                for (Double power: params.getPowerList()) {
                                    /*
                                     * add the sample size result to the list
                                     * if a failure occurs, an error code and message are
                                     * included with this object
                                     */
                                    results.add(getSampleSizeValue(params, test, method, alpha,
                                            sigmaScale, betaScale, power, quantile));
                                    if (Thread.currentThread().isInterrupted()) {
                                        return results;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return results;
    }

    /**
     * Runs a simulation to determine power values for the given
     * parameters.
     *
     * Note: quantile / unconditional power currently hard-coded to 1000
     * random instances of the baseline covariate x1000 error simulations.
     *
     * @see GLMMPowerParameters
     * @see GLMMPower
     * @param powerParams inputs to the simulation
     * @param iterations number of iterations to perform in the simulation
     * @return list of calculated power values.
     */
    @Override
    public List<Power> getSimulatedPower(PowerParameters powerParams, int iterations)
            throws PowerException {
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);

        // precalculate any computationally expensive matrices/constants,
        // update the parameters as needed - used for random covariates
        initialize(params);

        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();

        // calculate the power for all variations of the study design
        for (Test test: params.getTestList()) {
            for (PowerMethod method: params.getPowerMethodList()) {
                List<Double> quantileList = method == PowerMethod.QUANTILE_POWER
                                                ? params.getQuantileList()
                                                : NO_QUANTILES;

                for (Double alpha: params.getAlphaList()) {
                    for (Double sigmaScale: params.getSigmaScaleList()) {
                        for (Double betaScale: params.getBetaScaleList()) {
                            for (Double quantile: quantileList) {
                                for (Integer sampleSize: params.getSampleSizeList()) {
                                    /*
                                     * add the simulated power result to the list
                                     * if a failure occurs, an error code and message are
                                     * included with this object
                                     */
                                    results.add(getSimulatedPowerValue(params, test, method, alpha,
                                            sigmaScale, betaScale, sampleSize, quantile, iterations));
                                    if (Thread.currentThread().isInterrupted()) {
                                        return results;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return results;
    }

    /**
     * Find the best possible effect size (i.e. scale factor for the beta matrix)  to achieve
     * a specified power or list of powers.  Effect size is determined with a bisection search
     *
     * @see GLMMPowerParameters
     * @see GLMMPower
     * @param powerParams inputs to the effect size calculation
     * @return list of calculated power values.
     */
    @Override
    public List<Power> getDetectableDifference(PowerParameters powerParams)
            throws PowerException {
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);

        // precalculate any computationally expensive matrices/constants,
        // update the parameters as needed - used for random covariates
        initialize(params);

        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();

        // calculate the power for all variations of the study design
        for (Test test: params.getTestList()) {
            for (PowerMethod method: params.getPowerMethodList()) {
                List<Double> quantileList = method == PowerMethod.QUANTILE_POWER
                                                ? params.getQuantileList()
                                                : NO_QUANTILES;

                for (Double alpha: params.getAlphaList()) {
                    for (Double sigmaScale: params.getSigmaScaleList()) {
                        for (Integer sampleSize: params.getSampleSizeList()) {
                            for (Double quantile: quantileList) {
                                for (Double power: params.getPowerList()) {
                                    /*
                                     * add the detectable difference result to the list
                                     * if a failure occurs, an error code and message are
                                     * included with this object
                                     */
                                    results.add(getDetectableDifferenceValue(params,
                                            test, method, alpha, sigmaScale, power, sampleSize, quantile));
                                    if (Thread.currentThread().isInterrupted()) {
                                        return results;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return results;
    }

    /**
     * Perform any preliminary calculations / updates on the input matrices
     * @param params input parameters
     */
    private void initialize(GLMMPowerParameters params)
    {
        debug("entering initialize");

        // TODO: why isn't this done in PowerResourceHelper.studyDesignToPowerParameters?
        // if no power methods are specified, add conditional as a default
        if (params.getPowerMethodList().size() <= 0)
            params.addPowerMethod(PowerMethod.CONDITIONAL_POWER);

        // update sigma error and beta if we have a baseline covariate
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        debug("sigmaG", sigmaG);

        int numRandom = sigmaG != null ? sigmaG.getRowDimension() : 0;
        if (numRandom == 1)
        {
            // set the sigma error matrix to [sigmaY - sigmaYG * sigmaG-1 * sigmaGY]
            RealMatrix sigmaY = params.getSigmaOutcome();
            debug("sigmaY", sigmaY);

            RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
            debug("sigmaYG", sigmaYG);

            RealMatrix sigmaGY = sigmaYG.transpose();
            debug("sigmaYG transpose", sigmaGY);

            RealMatrix sigmaGInverse = new LUDecomposition(sigmaG).getSolver().getInverse();
            debug("sigmaG inverse", sigmaGInverse);

            RealMatrix sigmaError = forceSymmetric(sigmaY.subtract(sigmaYG.multiply(sigmaGInverse.multiply(sigmaGY))));
            debug("sigmaError = sigmaY - sigmaYG * sigmaG inverse * sigmaYG transpose", sigmaError);

            if (! MatrixUtils.isPositiveSemidefinite(sigmaError)) {
                throw new IllegalArgumentException(SIGMA_ERROR_NOT_POSITIVE_SEMIDEFINITE_MESSAGE);
            }
            params.setSigmaError(sigmaError);

            // calculate the betaG matrix and fill in the placeholder row for the random predictor
            FixedRandomMatrix beta = params.getBeta();
            beta.updateRandomMatrix(sigmaGInverse.multiply(sigmaGY));
            debug("beta with random part updated to sigmaG inverse * sigmaYG transpose", beta.getCombinedMatrix());
        }

        debug("exiting initialize");
    }

    /**
     * Ensure that all required matrices are specified, and that conformance is correct
     * @param params GLMM input parameters
     * @throws IllegalArgumentException
     */
    protected void validateMatrices(GLMMPowerParameters params)
            throws PowerException
    {
        if (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE
                && params.getSampleSizeForEstimates() <= params.getDesignMatrixRankForEstimates()) {
            throw new PowerException(NEGATIVE_NU_EST_MESSAGE, PowerErrorEnum.POWER_CI_NEGATIVE_NU_EST);
        }

        // convenience variables
        RealMatrix beta = params.getBeta().getCombinedMatrix();
        RealMatrix theta0 = params.getTheta();
        RealMatrix XEssence = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        int numRandom = sigmaG != null ? sigmaG.getRowDimension() : 0;

        // only allow at most 1 random predictor
        // TODO: handle multiple random predictors
        if (numRandom > 1)
            throw new PowerException("Too many random predictors - at most 1 is allowed",
                    PowerErrorEnum.MAX_RANDOM_PREDICTORS_EXCEEDED);

        // make sure all required matrices have been specified
        // note, we don't check U (within subject contrast), since it may be null in univariate cases
        if (beta == null)
            throw new PowerException("No beta (regression coefficients) matrix specified",
                    PowerErrorEnum.MISSING_MATRIX_BETA);
        if (XEssence == null)
            throw new PowerException("No design essence matrix specified",
                    PowerErrorEnum.MISSING_MATRIX_DESIGN);
        if (C == null)
            throw new PowerException("No between subject contrast (C) matrix specified",
                    PowerErrorEnum.MISSING_MATRIX_C);
        if (theta0 == null)
            throw new PowerException("No theta_null (null hypothesis) matrix specified",
                    PowerErrorEnum.MISSING_MATRIX_THETA_NULL);
        // create a default U if not specified
        if (U == null)
        {
            U = org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(beta.getColumnDimension());
            params.setWithinSubjectContrast(U);
        }

        // different variance/covariance matrices are specified depending on the presence
        // of random covariate
        if (numRandom == 0)
        {
            if (sigmaE == null)
                throw new PowerException("No sigma (error) matrix specified",
                        PowerErrorEnum.MISSING_MATRIX_SIGMA_E);
            if (!sigmaE.isSquare())
                throw new PowerException("Sigma error matrix must be square",
                        PowerErrorEnum.MATRIX_NONSQUARE_SIGMA_E);
            if (U.getRowDimension() != sigmaE.getRowDimension())
                throw new PowerException(
                        "Within subject contrast matrix "
                      + "(" + U.getRowDimension() + " x " + U.getColumnDimension() + ")"
                      + " does not conform with "
                      + "sigma matrix "
                      + "(" + sigmaE.getRowDimension() + " x " + sigmaE.getColumnDimension() + ")",
                        PowerErrorEnum.MATRIX_CONFORMANCE_U_SIGMA_E
                      );
        }
        else if (numRandom == 1)
        {
            // covariate (results not published for Wilk's Lambda or Pillai-Bartlett
            for (Test test: params.getTestList())
            {
                if (test != Test.HOTELLING_LAWLEY_TRACE && test != Test.UNIREP && test != Test.UNIREP_BOX &&
                        test != Test.UNIREP_GEISSER_GREENHOUSE && test != Test.UNIREP_HUYNH_FELDT)
                    throw new PowerException("With a random covariate, " +
                            "only Hotelling-Lawley and Unirep test statistics are supported",
                            PowerErrorEnum.UNKNOWN_TEST_REQUESTED_RANDOM);
            }

            if (sigmaG == null)
                throw new PowerException("No variance/covariance matrix specified for gaussian predictors",
                        PowerErrorEnum.MISSING_MATRIX_SIGMA_G);
            if (sigmaY == null)
                throw new PowerException("No variance/covariance matrix specified for response variables",
                        PowerErrorEnum.MISSING_MATRIX_SIGMA_Y);
            if (sigmaYG == null)
                throw new PowerException("No outcome / gaussian predictor covariance matrix specified",
                        PowerErrorEnum.MISSING_MATRIX_SIGMA_YG);

            // check conformance
            if (U.getRowDimension() != sigmaY.getRowDimension())
                throw new PowerException(
                        "Within subject contrast matrix "
                      + "(" + U.getRowDimension() + " x " + U.getColumnDimension() + ")"
                      + " does not conform with "
                      + "sigma matrix "
                      + "(" + sigmaY.getRowDimension() + " x " + sigmaY.getColumnDimension() + ")",
                        PowerErrorEnum.MATRIX_CONFORMANCE_U_SIGMA_Y
                      );
            if (sigmaG.getRowDimension() != sigmaYG.getColumnDimension())
                throw new PowerException("Outcome / Gaussian predictor covariance " +
                        "matrix does not conform with variance matrix for the gaussian predictor",
                        PowerErrorEnum.MATRIX_CONFORMANCE_SIGMA_G_SIGMA_YG);
            if (!sigmaY.isSquare())
                throw new PowerException("Variance/covariance matrix for " +
                        "response variables must be square",
                        PowerErrorEnum.MATRIX_NONSQUARE_SIGMA_Y);
            if (!sigmaG.isSquare())
                throw new PowerException("Variance/covariance matrix " +
                        "for gaussian predictors must be square",
                        PowerErrorEnum.MATRIX_NONSQUARE_SIGMA_G);
        }

        // check dimensionality
        if (C.getColumnDimension() != beta.getRowDimension())
            throw new PowerException("Between subject contrast matrix does not conform with beta matrix",
                    PowerErrorEnum.MATRIX_CONFORMANCE_C_BETA);
        if (beta.getColumnDimension() != U.getRowDimension())
            throw new PowerException("Within subject contrast matrix does not conform with beta matrix",
                    PowerErrorEnum.MATRIX_CONFORMANCE_BETA_U);
        if ((XEssence.getColumnDimension() != beta.getRowDimension() && numRandom == 0) ||
                (XEssence.getColumnDimension() + 1 != beta.getRowDimension() && numRandom > 0))
            throw new PowerException("Design matrix does not conform with beta matrix",
                    PowerErrorEnum.MATRIX_CONFORMANCE_X_BETA);
        if (C.getRowDimension() > C.getColumnDimension())
            throw new PowerException("Number of rows in between subject " +
                    "contrast must be less than or equal to the number of columns",
                    PowerErrorEnum.MATRIX_DIMENSION_C_TOO_MANY_ROWS);
        if (U.getColumnDimension() > U.getRowDimension())
            throw new PowerException("Number of columns in within " +
                    "subject contrast must be less than or equal to the number of rows",
                    PowerErrorEnum.MATRIX_DIMENSION_U_TOO_MANY_COLUMNS);
        if (theta0.getRowDimension() != C.getRowDimension())
            throw new PowerException("Number of rows in theta null " +
                    "must equal number of rows in between subject contrast",
                    PowerErrorEnum.MATRIX_CONFORMANCE_C_THETA_NULL);
        if (theta0.getColumnDimension() != U.getColumnDimension())
            throw new PowerException("Number of columns in theta null " +
                    "must equal number of columns in within subject contrast",
                    PowerErrorEnum.MATRIX_CONFORMANCE_U_THETA_NULL);

        // check rank of the design matrix
        int rankX = new SingularValueDecomposition(XEssence).getRank();
        if (rankX != Math.min(XEssence.getColumnDimension(), XEssence.getRowDimension()))
            throw new PowerException("Design matrix is not full rank: it is of rank " + rankX,
                    PowerErrorEnum.MATRIX_RANK_DESIGN_LTFR);

        // make sure design matrix is symmetric and positive definite
        // TODO: how to check this?
    }

    /**
     * Compute conditional power.  Conditional power is conditional on
     * a single instance of the design matrix, and is most appropriate for
     * designs with only categorical  predictors
     * @param params GLMN input parameters
     * @return conditional power
     */
    private double getConditionalPower(GLMMTest glmmTest, double alpha)
    {
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);

        // calculate the non-centrality parameter for the specified test statistic
        // under the null hypothesis
        double nonCentralityParam = glmmTest.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE);

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));
    }

    private static final double LOG2 = Math.log(2);

    /**
     * Calculate power by integrating over all possible values of the
     * non-centrality parameter.  Best used for designs with a
     * baseline covariate
     *
     * @param params GLMM input parameters
     * @return unconditional power
     * @throws PowerException
     */
    private double getUnconditionalPower(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist, double alpha)
            throws PowerException {
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);

        // get the distribution of the noncentrality parameter
        double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);

        debug("Fcrit = " + Fcrit + ", ndf = " + ndf + ", ddf = " + ddf);

        double h1 = nonCentralityDist.getH1();
        double h0 = nonCentralityDist.getH0();

        debug("h0 = " + h0 + ", h1 = " + h1);

        try {
            if (h1 < h0) {
                throw new IllegalArgumentException("Integration bounds are " + h0 + " and " + h1 + ".");
            }

            double integralResult;

            numberOfEvaluations = BigInteger.ZERO;

            if (h1 > h0) {
                // integrate over all values of non-centrality parameter from h0 to h1
                SimpsonIntegrator integrator = new SimpsonIntegrator();
                UnconditionalPowerIntegrand integrand = new UnconditionalPowerIntegrand(nonCentralityDist, Fcrit, ndf, ddf);

                integralResult = integrator.integrate(MAX_EVALUATIONS, integrand, h0, h1);
            } else {
                integralResult = 0;
            }

            final double I = integralResult;

            debug(new Supplier<Object>() {
                @Override
                public Object get() {
                    return
                        "done with integration: "
                      + "result = " + I + ", "
                      + "number of evaluations = " + numberOfEvaluations + ", "
                      + "log = " + (int) (Math.log(numberOfEvaluations.longValue() - 3)/LOG2 + 0.5);
                }
            });

            // create a noncentral F dist with non-centrality of H1
            NonCentralFDistribution fdist = new NonCentralFDistribution(ndf, ddf, h1);

            return (1 - fdist.cdf(Fcrit) - 0.5*integralResult);
        } catch (RuntimeException e) {
            LOGGER.warn("exiting getUnconditionalPower abnormally", e);

            throw new PowerException(e.getMessage(),
                    PowerErrorEnum.INTEGRATION_OVER_DISTRIBUTION_NONCENTRALITY_PARAMETER_FAILED);
        }
    }

    /**
     * Calculate quantile power by determining a specified quantile
     * of the non-centrality distribution.
     *
     * @param params GLMM input parameters
     * @return quantile power
     */
    private double getQuantilePower(GLMMTest glmmTest,
            NonCentralityDistribution nonCentralityDist, double alpha, double quantile)
    {
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);

        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        double nonCentralityParam = nonCentralityDist.inverseCDF(quantile);

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));
    }

    /**
     *  Find the sample size which achieves the desired power(s)
     *  specified in the input parameters.  Uses a bisection search.
     *
     * @param params GLMM input parameters
     * @return sample size
     */
    private GLMMPower getSampleSizeValue(GLMMPowerParameters params,
            Test test, PowerMethod method, double alpha,
            double sigmaScale, double betaScale, double targetPower, double quantile) {
        if (method == PowerMethod.UNCONDITIONAL_POWER) {
            GLMMPower power = new GLMMPower(test, alpha, targetPower, -1, -1,
                    betaScale, sigmaScale, method, quantile, null);
            power.setErrorMessage(
                "We have temporarily disabled sample size calculations using the unconditional power method "
              + "while we work on computational efficiency. "
              + "There are two alternatives."
              + "<ol>"
              + "<li>"
              + "You may perform a power calculation for a given sample size, "
              + "and iterate until you find a sample size and unconditional power that work for your design."
              + "</li>"
              + "<li>"
              + "You may calculate sample size using quantile power "
              + "calculated at the median power (0.50th quantile) instead of unconditional power. "
              + "As noted in Glueck and Muller (2003) "
              + "(see <a href=\"http://samplesizeshop.org/education/related-publications/\">Related Publications</a>), "
              + "quantile power is a very good approximation for unconditional power."
              + "</li>"
              + "</ol>"
              + "Thank you for your patience while we repair this functionality."
            );
            power.setErrorCode(PowerErrorEnum.POWER_METHOD_UNKNOWN);
            return power;
        }

        try {
            // rescale the beta and sigma matrices and create a test
            RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
            RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
            GLMMTest glmmTest = GLMMTestFactory.createGLMMTestForPower(test,
                    params.getFApproximationMethod(test),
                    params.getUnivariateCdfMethod(test),
                    params.getUnivariateEpsilonMethod(test),
                    params.getDesignEssence(),
                    params.getXtXInverse(),
                    STARTING_SAMPLE_SIZE,
                    params.getDesignRank(),
                    params.getBetweenSubjectContrast(),
                    params.getWithinSubjectContrast(),
                    params.getTheta(),
                    scaledBeta, scaledSigmaError,
                    (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE ?
                            params.getSampleSizeForEstimates() - params.getDesignMatrixRankForEstimates(): 0));

            // calculate the noncentrality distribution
            NonCentralityDistribution nonCentralityDist = null;
            if (method != PowerMethod.CONDITIONAL_POWER) {
                nonCentralityDist = new NonCentralityDistribution(test,
                        params.getDesignEssence(),
                        params.getXtXInverse(),
                        STARTING_SAMPLE_SIZE,
                        params.getBetweenSubjectContrast(),
                        params.getWithinSubjectContrast(),
                        params.getTheta(),
                        scaledBeta, scaledSigmaError,
                        params.getSigmaGaussianRandom(),
                        params.isNonCentralityCDFExact());
            }

            // Calculate the maximum valid per group N. This avoids multiplication which exceeeds
            // Integer.MAX_VALUE. Moreover, by limiting to MAX_SAMPLE_SIZE, we avoid calculations
            // that take many seconds to complete.
            int designEssenceRows = params.getDesignEssence().getRowDimension();
            int maxPerGroupN = Math.min(Integer.MAX_VALUE/designEssenceRows, MAX_SAMPLE_SIZE);

            /*
             * find the upper bound on sample size.  That is,
             * find a sample size which produces power greater than or
             * equal to the desired power.
             * If an error occurs, we set an error code in the power result
             */
            SampleSizeBound upperBound = getSampleSizeUpperBound(glmmTest, nonCentralityDist,
                    method, targetPower, alpha, quantile, maxPerGroupN);
            if (upperBound.getError() != null) {
                GLMMPower power = new GLMMPower(test, alpha, targetPower,
                        upperBound.getActualPower(), upperBound.getSampleSize(),
                        betaScale, sigmaScale, method, quantile, null);
                switch (upperBound.getError()) {
                case MAX_SAMPLE_SIZE_EXCEEDED:
                    power.setErrorMessage(
                     // TODO: is the following expression correct?
                     // "The total sample size for this case would exceed " + (maxPerGroupN * designEssenceRows) + ". "
                        "The total sample size for this case would be unreasonably large. "
                      + "For performance reasons, we are not computing its exact value."
                    );
                    power.setErrorCode(PowerErrorEnum.MAX_SAMPLE_SIZE_EXCEEDED);
                    break;
                case SAMPLE_SIZE_UNDEFINED:
                    power.setErrorMessage(NO_MEAN_DIFFERENCE_MESSAGE);
                    power.setErrorCode(PowerErrorEnum.SAMPLE_SIZE_UNDEFINED);
                    break;
                case SAMPLE_SIZE_UNDEFINED_DUE_TO_EXCEPTION:
                    power.setErrorMessage("Sample size not well defined.");
                    power.setErrorCode(PowerErrorEnum.SAMPLE_SIZE_UNDEFINED);
                    break;
                }
                return power;
            }

            /*
             *  find the lower bound on sample size
             *  We search for a lower bound since the F approximations
             *  may be unstable for very small samples
             */
            // TODO: isn't the lower bound just half of the upper bound?????
            SampleSizeBound lowerBound = getSampleSizeLowerBound(glmmTest, nonCentralityDist,
                    method, upperBound, alpha, quantile);
            assert lowerBound.getError() == null;

            /*
             * At this point we have valid boundaries for searching.
             * There are two possible scenarios
             * 1. The upper bound == lower bound.
             * 2. The upper bound != lower bound and lower bound exceeds required power.
             * In this case we just take the value at the lower bound.
             * 3. The upper bound != lower bound and lower bound is less than the required power.
             * In this case we bisection search
             */
            double calculatedPower;
            int perGroupSampleSize;

            if (upperBound.getSampleSize() == lowerBound.getSampleSize()) {
                // case 1
                calculatedPower = upperBound.getActualPower();
                perGroupSampleSize = upperBound.getSampleSize();
            } else if (lowerBound.getActualPower() >= targetPower) {
                // case 2
                calculatedPower = lowerBound.getActualPower();
                perGroupSampleSize = lowerBound.getSampleSize();
            } else {
                // case 3, bisection search time!
                // create a bisection search function to find the best per group sample size
                BisectionSolver solver =  new BisectionSolver();
                SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(glmmTest, nonCentralityDist,
                        method, targetPower, alpha, quantile);
                double solution = solver.solve(MAX_ITERATIONS, sampleSizeFunc,
                        lowerBound.getSampleSize(), upperBound.getSampleSize());
                perGroupSampleSize = (int) Math.rint(solution); // see https://samplesizeshop.atlassian.net/browse/SSS-120
                glmmTest.setPerGroupSampleSize(perGroupSampleSize);
                if (nonCentralityDist != null) nonCentralityDist.setPerGroupSampleSize(perGroupSampleSize);
                calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
            }

            // build a confidence interval if requested
            GLMMPowerConfidenceInterval ci;
            if (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE) {
                ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
                        params.getAlphaLowerConfidenceLimit(),
                        params.getAlphaUpperConfidenceLimit(),
                        params.getSampleSizeForEstimates(),
                        params.getDesignMatrixRankForEstimates(),
                        alpha, glmmTest);
            } else {
                ci = null;
            }

            return new GLMMPower(test, alpha, targetPower, calculatedPower,
                    MatrixUtils.getTotalSampleSize(params.getDesignEssence(), perGroupSampleSize),
                    betaScale, sigmaScale, method, quantile, ci);
        } catch (PowerException pe) {
            GLMMPower powerValue = new GLMMPower(test, alpha, targetPower, -1, -1,
                    betaScale, sigmaScale, method, quantile, null);
            powerValue.setErrorCode(pe.getErrorCode());
            powerValue.setErrorMessage(pe.getMessage());
            return powerValue;
        }
    }

    /**
     *
     * @param glmmTest the statistical test
     * @param nonCentralityDist the noncentrality distribution if
     * the design includes a covariate
     * @param method power calculation method
     * @param alpha
     * @param quantile
     * @return
     */
    private SampleSizeBound getSampleSizeLowerBound(GLMMTest glmmTest,
            NonCentralityDistribution nonCentralityDist, PowerMethod method,
            SampleSizeBound upperBound, double alpha, double quantile) {
        assert upperBound.getError() == null;

        int upperN = upperBound.getSampleSize();
        if (upperN <= 2) {
            /* If the upper bound is 2, then there is no smaller sample size
             * to use as a lower bound
             *  If the upper bound is -1, then we could not find a valid sample size
             *  which exceeds the desired power
             *  In either case, we just return a lower bound equivalent to the
             *  upper bound
             */
            return new SampleSizeBound(upperN, upperBound.getActualPower());
        } else {
            /* we have a valid upper bound, so search for the smallest valid sample
             * size for this design */
            double calculatedPower = 0;
            // this becomes 2 as we enter the loop, since N=1 makes no sense statistically
            int lowerBound = 1;
            do {
                lowerBound++;
                try {
                    glmmTest.setPerGroupSampleSize(lowerBound);
                    if (nonCentralityDist != null) {
                        nonCentralityDist.setPerGroupSampleSize(lowerBound);
                    }
                    calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
                    // if we don't throw an exception, then we have a valid minimum
                    break;
                } catch (Exception e) {
                    // just keep iterating until we find a minimum valid sample size
                    LOGGER.info("Ignoring exception as we iterate to find a valid minimum:", e);
                }
            } while (lowerBound < upperN && !Thread.currentThread().isInterrupted());

            return new SampleSizeBound(lowerBound, calculatedPower);
        }
    }

    /**
     * Returns true if all values in the beta matrix are identical
     * @param beta beta matrix
     * @return true if no mean difference, false otherwise
     */
    private boolean noMeanDifference(GLMMTest test) {
        // get the difference between theta null and the alternative
        RealMatrix sumSqHypothesis = test.getHypothesisSumOfSquares();
        // check if there is at least one non-zero value
        if (sumSqHypothesis != null) {
            for(int r = 0; r < sumSqHypothesis.getRowDimension(); r++) {
                for(int c = 0; c < sumSqHypothesis.getColumnDimension(); c++) {
                    if (Math.abs(sumSqHypothesis.getEntry(r, c)) > EPSILON) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Determine the upper bound for the bisection search used in
     * calculation of sample size
     *
     * @param params GLMM input parameters
     * @return upper bound on sample size to achieve desired power
     */
    private SampleSizeBound getSampleSizeUpperBound(GLMMTest glmmTest,
            NonCentralityDistribution nonCentralityDist,
            PowerMethod method, double targetPower,
            double alpha, double quantile, int maxPerGroupN) {
        // check if no mean difference.  In this case, sample size is undefined and
        // power is always alpha
        if (noMeanDifference(glmmTest)) {
            return new SampleSizeBound(-1, alpha, SampleSizeError.SAMPLE_SIZE_UNDEFINED);
        }

        // otherwise, keep ramping up sample size until we exceed the desired power
        int upperBound = STARTING_SAMPLE_SIZE;
        double currentPower;
        do {
            upperBound += upperBound;
            // check for overflows
            if (upperBound > maxPerGroupN) {
                upperBound = maxPerGroupN;
            }
            try {
                glmmTest.setPerGroupSampleSize(upperBound);
                if (nonCentralityDist != null) {
                    nonCentralityDist.setPerGroupSampleSize(upperBound);
                }
                currentPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
            } catch (Exception e) {
                // ignore steps which yield invalid degrees of freedom
                LOGGER.warn("exception getting power by type", e);
                return new SampleSizeBound(-1, alpha, SampleSizeError.SAMPLE_SIZE_UNDEFINED_DUE_TO_EXCEPTION);
            }
        } while (currentPower <= targetPower && upperBound < maxPerGroupN && !Thread.currentThread().isInterrupted());
        if (currentPower < targetPower) {
            // no sample size meets the criteria, so return an error
            return new SampleSizeBound(-1, alpha, SampleSizeError.MAX_SAMPLE_SIZE_EXCEEDED);
        } else {
            return new SampleSizeBound(upperBound, currentPower);
        }
    }

    /**
     * Perform a bisection search to determine effect size
     * @param params GLMM input parameters
     * @return detectable difference object
     * @throws IllegalArgumentException
     * @throws ArithmeticException
     */
    private GLMMPower getDetectableDifferenceValue(GLMMPowerParameters params,
            Test test, PowerMethod method, double alpha,
            double sigmaScale, double targetPower, int sampleSize, double quantile)
                    throws PowerException {
        RealMatrix scaledBeta = params.getBeta().scalarMultiply(STARTING_BETA_SCALE, true);
        RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTestForPower(test,
                params.getFApproximationMethod(test),
                params.getUnivariateCdfMethod(test),
                params.getUnivariateEpsilonMethod(test),
                params.getDesignEssence(),
                params.getXtXInverse(),
                sampleSize,
                params.getDesignRank(),
                params.getBetweenSubjectContrast(),
                params.getWithinSubjectContrast(),
                params.getTheta(),
                scaledBeta, scaledSigmaError,
                (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE ?
                        params.getSampleSizeForEstimates() - params.getDesignMatrixRankForEstimates(): 0));

        NonCentralityDistribution nonCentralityDist = null;
        if (method != PowerMethod.CONDITIONAL_POWER)
        {
            nonCentralityDist = new NonCentralityDistribution(test,
                    params.getDesignEssence(),
                    params.getXtXInverse(),
                    sampleSize,
                    params.getBetweenSubjectContrast(),
                    params.getWithinSubjectContrast(),
                    params.getTheta(),
                    scaledBeta, scaledSigmaError,
                    params.getSigmaGaussianRandom(),
                    params.isNonCentralityCDFExact());
        }

        BisectionSolver solver =  new BisectionSolver();

        DetectableDifferenceFunction ddFunc = new DetectableDifferenceFunction(glmmTest,
                nonCentralityDist, method, targetPower, alpha, quantile, params.getBeta());

        double lowerBound = 1.0E-10; // need non-zero lower bound or noncentrality dist malfunctions

        // check if the lower bound beta scale is already greater than the target power
        scaledBeta = params.getBeta().scalarMultiply(lowerBound, true);
        glmmTest.setBeta(scaledBeta);
        if (nonCentralityDist != null) nonCentralityDist.setBeta(scaledBeta);
        double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
        // find the detectable difference (i.e. beta scale) if the lower bound does not exceed the target power
        double betaScale = lowerBound;
        if (calculatedPower < targetPower)
        {
            DetectableDifferenceBound upperBound =
                    getDetectableDifferenceUpperBound(glmmTest,
                            nonCentralityDist, params.getBeta(),method,
                            targetPower, alpha, quantile);
            if (upperBound.getError() != null) {
                GLMMPower power = new GLMMPower(test, alpha, targetPower,
                        upperBound.getActualPower(),
                        MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                        upperBound.getBetaScale(),
                        sigmaScale, method, quantile, null);
                switch (upperBound.getError()) {
                case MAX_BETA_SCALE_EXCEEDED:
                    power.setErrorMessage("Failed to find valid upper bound on beta scale.");
                    power.setErrorCode(PowerErrorEnum.MAX_BETA_SCALE_EXCEEDED);
                    break;
                case BETA_SCALE_UNDEFINED:
                    power.setErrorMessage("Beta scale not well defined for no difference.");
                    power.setErrorCode(PowerErrorEnum.BETA_SCALE_UNDEFINED);
                    break;
                }
                return power;
            }
            betaScale = solver.solve(MAX_ITERATIONS, ddFunc, lowerBound,
                    upperBound.getBetaScale());
        }
        // calculate actual power associated with this beta scale
        scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
        glmmTest.setBeta(scaledBeta);
        if (nonCentralityDist != null) nonCentralityDist.setBeta(scaledBeta);
        calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);

        // build a confidence interval if requested
        GLMMPowerConfidenceInterval ci;
        if (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE) {
            ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
                    params.getAlphaLowerConfidenceLimit(),
                    params.getAlphaUpperConfidenceLimit(),
                    params.getSampleSizeForEstimates(),
                    params.getDesignMatrixRankForEstimates(),
                    alpha, glmmTest);
        } else {
            ci = null;
        }

        return new GLMMPower(test, alpha, targetPower, calculatedPower,
                MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                betaScale, sigmaScale, method, quantile, ci);
    }

    /**
     * Get the upper bound for the bisection search used to determine effect size
     * @param params GLMM input parameters
     * @return upper bound on beta scale
     */
    private DetectableDifferenceBound getDetectableDifferenceUpperBound(GLMMTest glmmTest,
            NonCentralityDistribution nonCentralityDist, FixedRandomMatrix beta,
            PowerMethod method, double targetPower, double alpha, double quantile)
    {
        double maxBetaScale =
                Double.MAX_VALUE / MatrixUtils.getMaxValue(beta.getCombinedMatrix());
        double upperBound = STARTING_BETA_SCALE;
        double prevBound;
        double currentPower = 0.0;
        do {
            prevBound = upperBound;
            upperBound *= 2;
            // check for double overflows
            if (upperBound < prevBound || upperBound > maxBetaScale) {
                upperBound = maxBetaScale;
            }
            try {
                RealMatrix scaledBeta = beta.scalarMultiply(upperBound, true);
                glmmTest.setBeta(scaledBeta);
                if (nonCentralityDist != null) nonCentralityDist.setBeta(scaledBeta);
                currentPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
            } catch (Exception e) {
                // ignore steps which yield invalid degrees of freedom
            }
        } while (currentPower <= targetPower && upperBound <= Double.MAX_VALUE && !Thread.currentThread().isInterrupted());
        if (currentPower < targetPower) {
            // no sample size meets the criteria, so return an error
            return new DetectableDifferenceBound(Double.NaN, alpha,
                    DetectableDifferenceError.MAX_BETA_SCALE_EXCEEDED);
        } else {
            return new DetectableDifferenceBound(upperBound, currentPower);
        }
    }


    //    /**
    //     * Simulate power for the general linear multivariate model based on
    //     * the input matrices.
    //     *
    //     * @param params Container for input matrices
    //     * @param iterations number of simulation samples/iterations
    //     * @return simulated power value
    //     */
    //    public double simulatePower(GLMMPowerParameters params, int iterations)
    //            throws IllegalArgumentException
    //    {
    //        // get total observations, N, and rank of design matrix
    //        RealMatrix XEssence = params.getDesignEssence();
    //        double N = XEssence.getRowDimension()*params.getCurrentSampleSize();
    //        double rankX = params.getDesignRank();
    //
    //        // create a normal distribution for generating random errors
    //        Normal normalDist = new Normal();
    //        normalDist.setSeed(1234);
    //        int rejectionCount = 0;
    //        // create an error matrix here, so we don't have to reallocate every time
    //        Array2DRowRealMatrix randomNormals =
    //            new Array2DRowRealMatrix((int) N, params.getScaledBeta().getColumnDimension());
    //
    //        if (params.getSigmaGaussianRandom() != null &&
    //                params.getCurrentPowerMethod() != PowerMethod.CONDITIONAL_POWER)
    //        {
    //            double[] powerValues = new double[SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL];
    //
    //            for(int gInstance = 0; gInstance < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; gInstance++)
    //            {
    //                // force a new realization of the design matrix (i.e. a new covariate column)
    //                RealMatrix X = getFullDesignMatrix(params.getDesignEssence(),
    //                        params.getSigmaGaussianRandom(),    perGroupSize);
    //                rejectionCount = 0;
    //                for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
    //                {
    //                    GLMMPowerParameters params, Test test,
    //                    RealMatrix X, RealMatrix XtXinverse,
    //                    RealMatrix scaledBeta,
    //                    RandomErrorMatrix randomErrors,
    //                    int perGroupSampleSize, double alpha
    //
    //
    //                    if (simulateAndFitModel(params, test, X, null, scaledBeta, random normalDist, randomNormals, N, rankX)) rejectionCount++;
    //                }
    //                powerValues[gInstance] = (((double) rejectionCount) /
    //                        ((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
    //            }
    //
    //            switch (params.getCurrentPowerMethod())
    //            {
    //            case UNCONDITIONAL_POWER:
    //                return StatUtils.mean(powerValues);
    //            case QUANTILE_POWER:
    //                return StatUtils.percentile(powerValues, params.getCurrentQuantile());
    //            default:
    //                throw new IllegalArgumentException("Unknown power method");
    //            }
    //        }
    //        else
    //        {
    //            // run the simulations
    //            for(int i = 0; i < iterations; i++)
    //            {
    //                if (simulateAndFitModel(params, normalDist, randomNormals, N, rankX)) rejectionCount++;
    //            }
    //            return ((double) rejectionCount) / ((double) iterations);
    //        }
    //    }

    /**
     * Simulate the error matrix to generate a single realization of the data, then
     * fit the model and determine if the null hypothesis is rejected
     *
     * @param params GLMM parameter set
     * @param normalDist normal distribution object for generating random errors
     * @param N total number of observations
     * @param rankX rank of the design matrix
     * @return true if null is rejected, false otherwise
     */
    private SimulationFit simulateAndFitModel(GLMMPowerParameters params, Test test,
            RealMatrix X, RealMatrix XtXinverse,
            RealMatrix scaledBeta,
            RandomErrorMatrix randomErrors,
            int perGroupSampleSize, double alpha)
    {
        // get a new set of errors
        RealMatrix error = randomErrors.random();

        // calculate simulated Y based on Y = X beta + error
        RealMatrix Ysim = (X.multiply(scaledBeta)).add(error);

        // build a test object for the simulated matrices
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTestForDataAnalysis(test,
                params.getFApproximationMethod(test),
                params.getUnivariateCdfMethod(test),
                params.getUnivariateEpsilonMethod(test),
                X, XtXinverse, params.getDesignRank(), Ysim,
                params.getBetweenSubjectContrast().getCombinedMatrix(),
                params.getWithinSubjectContrast(),
                params.getTheta());
        ModelFit fit = glmmTest.getModelFit();

        // return the fit information for the simulated matrices
        return new SimulationFit(fit.Pvalue, fit.Fvalue, fit.numeratorDF,
                fit.denominatorDF, Ysim, X.multiply(fit.beta), fit.sigma, fit.beta);
    }

    /**
     * Returns a simulated power sample for each combination of parameters for a
     * design with a random covariate (GLMM(F,g)).  Currently produces
     * a sample of size 1000
     *
     * @param params power parameters
     * @param size size of the sample
     * @return list of simulated power samples
     * @throws IllegalArgumentException
     */
    public List<SimulatedPower[]> getSimulatedPowerSample(GLMMPowerParameters params, int size)
            throws PowerException {
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);

        // precalculate any computationally expensive matrices/constants,
        // update the parameters as needed - used for random covariates
        initialize(params);

        if (params == null || params.getSigmaGaussianRandom() == null)
            throw new IllegalArgumentException("Power samples can only be generated for designs with random predictors");
        if (size <= 0) throw new IllegalArgumentException("Iterations must be positive");

        // list of simulated power results
        ArrayList<SimulatedPower[]> results = new ArrayList<SimulatedPower[]>();

        // calculate the simulated power for all variations of the study design
        for (Test test: params.getTestList()) {
            for (PowerMethod method: params.getPowerMethodList()) {
                // only generate samples for quantile or unconditional power
                if (method == PowerMethod.CONDITIONAL_POWER) continue;

                List<Double> quantileList = method == PowerMethod.QUANTILE_POWER
                                                ? params.getQuantileList()
                                                : NO_QUANTILES;

                for (Double alpha: params.getAlphaList()) {
                    for (Double sigmaScale: params.getSigmaScaleList()) {
                        for (Double betaScale: params.getBetaScaleList()) {
                            for (Double quantile: quantileList) {
                                for (Integer sampleSize: params.getSampleSizeList()) {
                                    /*
                                     * add the simulated power sample result to the list
                                     * if a failure occurs, an error code and message are
                                     * included with this object
                                     */
                                    results.add(getSimulatedPowerSampleValue(params, test, method, alpha,
                                            sigmaScale, betaScale, sampleSize, quantile, size));
                                    if (Thread.currentThread().isInterrupted()) {
                                        return results;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return results;
    }

    /**
     * Generate a sample of powers for a design with a random covariate
     * @param params power parameters
     * @param iterations size of the sample
     * @return
     */
    private SimulatedPower[] getSimulatedPowerSampleValue(GLMMPowerParameters params,
            Test test, PowerMethod method, double alpha,
            double sigmaScale, double betaScale, int sampleSize, double quantile, int iterations)
    {
        // scale beta and sigma
        RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
        RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
        // get random errors
        RandomErrorMatrix randomErrors =
                new RandomErrorMatrix(MatrixUtils.getTotalSampleSize(params.getDesignEssence(),
                        sampleSize), scaledBeta.getColumnDimension(), scaledSigmaError);
        randomErrors.setPositivityThreshold(positivityThreshold);
        randomErrors.setSymmetryThreshold(symmetryThreshold);
        SimulatedPower[] simPower = new SimulatedPower[iterations];

        for(int gInstance = 0; gInstance < iterations; gInstance++)
        {
            // force a new realization of the design matrix (i.e. a new covariate column)
            RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(),
                    params.getSigmaGaussianRandom(), sampleSize, seed);
            RealMatrix XtXinverse = new LUDecomposition(X.transpose().multiply(X)).getSolver().getInverse();

            int rejectionCount = 0;
            simPower[gInstance] =
                    new SimulatedPower(scaledBeta.getRowDimension(), scaledBeta.getColumnDimension());
            for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
            {
                SimulationFit fit =
                        simulateAndFitModel(params, test, X, XtXinverse, scaledBeta, randomErrors,
                                sampleSize, alpha);
                if (fit.Pvalue <= alpha) rejectionCount++;
                simPower[gInstance].accumulateBeta(fit.beta);
            }
            simPower[gInstance].setPower(((double) rejectionCount) /
                    ((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
        }

        return simPower;
    }

    /**
     * Convenience function to determine which power method to use
     *
     * @param params GLMM input parameters
     * @return conditional/quantile/unconditional power
     */
    private GLMMPower getPowerValue(GLMMPowerParameters params,
            Test test, PowerMethod method, double alpha,
            double sigmaScale, double betaScale, int sampleSize, double quantile) {
        if (method == PowerMethod.UNCONDITIONAL_POWER && sampleSize > MAX_SAMPLE_SIZE_FOR_UNCONDITIONAL_POWER) {
            GLMMPower power = new GLMMPower(test, alpha, -1, -1,
                    MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                    betaScale, sigmaScale, method, quantile, null);
            power.setErrorMessage(
                "For power calculations using the unconditional power method, we require that the Smallest "
              + "Group Size not exceed " + MAX_SAMPLE_SIZE_FOR_UNCONDITIONAL_POWER + ", for "
              + "performance reasons. "
              + "If your Smallest Group Size does exceed " + MAX_SAMPLE_SIZE_FOR_UNCONDITIONAL_POWER + ", "
              + "you may instead calculate power using quantile power "
              + "calculated at the median power (0.50th quantile) instead of unconditional power. "
              + "As noted in Glueck and Muller (2003) "
              + "(see <a href=\"http://samplesizeshop.org/education/related-publications/\">Related Publications</a>), "
              + "quantile power is a very good approximation for unconditional power."
            );
            power.setErrorCode(PowerErrorEnum.POWER_METHOD_UNKNOWN);
            return power;
        }

        try {
            RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
            RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
            GLMMTest glmmTest = GLMMTestFactory.createGLMMTestForPower(test,
                    params.getFApproximationMethod(test),
                    params.getUnivariateCdfMethod(test),
                    params.getUnivariateEpsilonMethod(test),
                    params.getDesignEssence(),
                    params.getXtXInverse(),
                    sampleSize,
                    params.getDesignRank(),
                    params.getBetweenSubjectContrast(),
                    params.getWithinSubjectContrast(),
                    params.getTheta(),
                    scaledBeta, scaledSigmaError,
                    (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE ?
                            params.getSampleSizeForEstimates() - params.getDesignMatrixRankForEstimates(): 0));

            // check if no mean difference.  In this case,
            // power is always alpha
            if (noMeanDifference(glmmTest)) {
                GLMMPower power = new GLMMPower(test, alpha, alpha, alpha,
                                   MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                                   betaScale, sigmaScale, method, quantile, null);
                power.setErrorMessage(NO_MEAN_DIFFERENCE_MESSAGE);
                power.setErrorCode(PowerErrorEnum.SAMPLE_SIZE_UNDEFINED);
                return power;
            }

            NonCentralityDistribution nonCentralityDist = null;
            if (method != PowerMethod.CONDITIONAL_POWER)
            {
                nonCentralityDist = new NonCentralityDistribution(
                        test,
                        params.getDesignEssence(),
                        params.getXtXInverse(),
                        sampleSize,
                        params.getBetweenSubjectContrast(),
                        params.getWithinSubjectContrast(),
                        params.getTheta(),
                        scaledBeta,
                        scaledSigmaError,
                        params.getSigmaGaussianRandom(),
                        params.isNonCentralityCDFExact());
            }

            // calculate the power
            double power = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);

            // build a confidence interval if requested
            GLMMPowerConfidenceInterval ci;
            if (params.getConfidenceIntervalType() != ConfidenceIntervalType.NONE) {
                ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
                        params.getAlphaLowerConfidenceLimit(),
                        params.getAlphaUpperConfidenceLimit(),
                        params.getSampleSizeForEstimates(),
                        params.getDesignMatrixRankForEstimates(),
                        alpha, glmmTest);
            } else {
                ci = null;
            }

            return new GLMMPower(test, alpha, power, power,
                    MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                    betaScale, sigmaScale, method, quantile, ci);
        } catch (PowerException pe) {
            GLMMPower powerValue = new GLMMPower(test, alpha, -1, -1,
                    MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
                    betaScale, sigmaScale, method, quantile, null);
            powerValue.setErrorCode(pe.getErrorCode());
            powerValue.setErrorMessage(pe.getMessage());
            return powerValue;
        }
    }

    private GLMMPower getSimulatedPowerValue(GLMMPowerParameters params,
            Test test, PowerMethod method, double alpha,
            double sigmaScale, double betaScale, int sampleSize, double quantile,
            int iterations)
    {
        // scale beta and sigma
        RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
        RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
        // get random errors
        RandomErrorMatrix randomErrors =
                new RandomErrorMatrix(MatrixUtils.getTotalSampleSize(params.getDesignEssence(),
                        sampleSize), scaledSigmaError.getColumnDimension(), scaledSigmaError);
        randomErrors.setPositivityThreshold(positivityThreshold);
        randomErrors.setSymmetryThreshold(symmetryThreshold);

        double power = Double.NaN;

        if (params.getSigmaGaussianRandom() != null &&
                method != PowerMethod.CONDITIONAL_POWER)
        {
            double[] powerValues = new double[SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL];

            for(int gInstance = 0; gInstance < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; gInstance++)
            {
                int rejectionCount = 0;
                // force a new realization of the design matrix (i.e. a new covariate column)
                RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(),
                        params.getSigmaGaussianRandom(), sampleSize, seed);
                RealMatrix XtXinverse = new LUDecomposition(X.transpose().multiply(X)).getSolver().getInverse();
                for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
                {
                    SimulationFit fit =
                            simulateAndFitModel(params, test, X, XtXinverse, scaledBeta, randomErrors,
                                    sampleSize, alpha);
                    if (fit.Pvalue <= alpha) rejectionCount++;
                }
                powerValues[gInstance] = (((double) rejectionCount) /
                        ((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
            }

            switch (method)
            {
            case UNCONDITIONAL_POWER:
                power = StatUtils.mean(powerValues);
                break;
            case QUANTILE_POWER:
                power= StatUtils.percentile(powerValues, quantile);
                break;
            default:
                throw new IllegalArgumentException("Unknown power method");
            }
        }
        else
        {
            // GLMM(F) design, conditional power
            // run the simulations
            RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(), sampleSize);
            RealMatrix XtXInverse = params.getXtXInverse().scalarMultiply(1/(double)sampleSize);
            int rejectionCount = 0;
            for(int i = 0; i < iterations; i++)
            {
                SimulationFit fit =
                        simulateAndFitModel(params, test, X, XtXInverse, scaledBeta, randomErrors,
                                sampleSize, alpha);
                if (fit.Pvalue <= alpha) rejectionCount++;
            }
            power =  ((double) rejectionCount) / ((double) iterations);
        }

        return new GLMMPower(test, alpha, power, power,
                MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize), betaScale,
                sigmaScale, method, quantile, null);
    }

    /**
     * Get an individual power instance
     * @param glmmTest
     * @param nonCentralityDist
     * @param method
     * @param targetPower
     * @param alpha
     * @param quantile
     * @return
     */
    private double getPowerByType(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
                                    PowerMethod method, double alpha, double quantile)
            throws PowerException {
        // calculate the power
        double power;

        switch (method) {
        case QUANTILE_POWER:
            power = getQuantilePower(glmmTest, nonCentralityDist, alpha, quantile);
            break;

        case UNCONDITIONAL_POWER:
            power = getUnconditionalPower(glmmTest, nonCentralityDist, alpha);
            break;

        case CONDITIONAL_POWER:
        default:
            power = getConditionalPower(glmmTest, alpha);
            break;
        }

        return power;
    }

    /**
     * Set the positivity threshold for Cholesky decomposition during simulation.  This
     * allows Cholesky decomposition for matrices with very small negative values.
     */
    public void setPositivityThreshold(double threshold)
    {
        this.positivityThreshold = threshold;
    }

    /**
     * Set the symmetry threshold for Cholesky decomposition during simulation.  This
     * allows Cholesky decomposition for matrices with very small differences between
     * symmetric cells.
     */
    public void setSymmetryThreshold(double threshold)
    {
        this.symmetryThreshold = threshold;
    }

    /**
     * A convenience method for DEBUG logging of a message.
     *
     * @param message The message.
     */
    private static void debug(Object message) {
        LOGGER.debug(message);
    }

    /**
     * A convenience method for DEBUG logging of a supplied message.
     *
     * @param supplier The message supplier.
     */
    private static void debug(Supplier<Object> supplier) {
        LOGGER.debug(supplier);
    }

    /**
     * A convenience method for DEBUG logging of a matrix
     * with a label.
     *
     * @param label      The label.
     * @param realMatrix The matrix.
     */
    private static void debug(String label, RealMatrix realMatrix) {
        LOGGER.debug(MatrixUtilities.logMessageSupplier(label, realMatrix));
    }
}
