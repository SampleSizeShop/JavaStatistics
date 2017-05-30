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
package edu.cudenver.bios.power.glmm;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.distribution.ChiSquareTerm;
import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.distribution.WeightedSumOfNoncentralChiSquaresDistribution;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtilities;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.power.PowerErrorEnum;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.utils.Logger;

import static edu.cudenver.bios.matrix.MatrixUtilities.forceSymmetric;

/**
 * Class representing the distribution of the non-centrality parameter in
 * the general linear multivariate model.  Used by the GLMMPowerCalculator class
 * for computing unconditional and quantile power.
 *
 * @see edu.cudenver.bios.power.GLMMPowerCalculator
 * @author Sarah Kreidler
 */
public class NonCentralityDistribution
{
    private static final String NOT_POSITIVE_DEFINITE =
            "Unfortunately, there is no solution for this combination of input parameters. "
        +   "A matrix that arose during the computation is not positive definite. "
        +   "It may be possible to reduce expected covariate/response correlations "
        +   "and obtain a soluble combination."
        ;

    private static final int MAX_ITERATIONS = 10000;
    private static final double ACCURACY = 0.001;
    // intermediate forms
    protected RealMatrix T1 = null;
    protected RealMatrix FT1 = null;
    protected RealMatrix S = null;
    protected RealMatrix mzSq = null;
    protected double H1;
    protected double H0 = 0;
    int qF;
    int a;
    double N;
    double[] sEigenValues;
    int sStar = 0;
    // indicates if an "exact" cdf should be calculated via Davie's algorithm or
    // with the Satterthwaite approximation from Glueck & Muller
    protected boolean exact;

    // cache input parameters - needed for dynamic reset of sample size and beta matrix
    protected Test test;
    protected RealMatrix FEssence;
    protected RealMatrix FtFinverse;
    protected int perGroupN;
    protected FixedRandomMatrix CFixedRand;
    protected RealMatrix U;
    protected RealMatrix thetaNull;
    protected RealMatrix beta;
    protected RealMatrix sigmaError;
    protected RealMatrix sigmaG;

    private static final Logger LOGGER = Logger.getLogger(NonCentralityDistribution.class);

    /**
     * Function calculating the difference between the probability of a target quantile
     * and the  (used by the bisection solver from Apache Commons Math)
     * @see org.apache.commons.math.analysis.UnivariateRealFunction
     */
    private class NonCentralityQuantileFunction implements UnivariateFunction
    {
        protected double quantile;

        public NonCentralityQuantileFunction(double quantile)
        {
            this.quantile = quantile;
        }

        public double value(double n)
        {
            try {
                return cdf(n) - quantile;
            } catch (PowerException pe) {
                throw new IllegalArgumentException(pe.getMessage(), pe);
            }
        }
    }

    /**
     * Create a non-centrality distribution for the specified inputs.
     * @param params GLMM input parameters
     * @param exact if true, Davie's algorithm will be used to compute the cdf,
     * otherwise a Satterthwaite style approximation is used.
     * @throws IllegalArgumentException
     */
    public NonCentralityDistribution(Test test, RealMatrix FEssence, RealMatrix FtFinverse, int perGroupN,
            FixedRandomMatrix CFixedRand, RealMatrix U,
            RealMatrix thetaNull, RealMatrix beta,
            RealMatrix sigmaError, RealMatrix sigmaG, boolean exact)
    throws PowerException
    {
        debug("CREATING NonCentralityDistribution");
        debug("begin parameters");
        debug("test: " + test);
        debug("FEssence:", FEssence);
        debug("FtFinverse:", FtFinverse);
        debug("perGroupN: " + perGroupN);
        debug("CFixedRand:", CFixedRand.getCombinedMatrix());
        debug("U:", U);
        debug("thetaNull:", thetaNull);
        debug("beta:", beta);
        debug("sigmaError:", sigmaError);
        debug("sigmaG:", sigmaG);
        debug("exact: " + exact);
        debug("end parameters");

        initialize(test, FEssence, FtFinverse, perGroupN, CFixedRand, U, thetaNull, beta,
                sigmaError, sigmaG, exact);
    }

    /**
     * Pre-calculate intermediate matrices, perform setup, etc.
     */
    private void initialize(Test test, RealMatrix FEssence, RealMatrix FtFinverse, int perGroupN,
            FixedRandomMatrix CFixedRand, RealMatrix U,
            RealMatrix thetaNull, RealMatrix beta,
            RealMatrix sigmaError, RealMatrix sigmaG, boolean exact)
    throws PowerException {
        debug("entering initialize");

        // reset member variables
        this.T1 = null;
        this.FT1 = null;
        this.S = null;
        this.mzSq = null;
        this.H0 = 0;
        this.sStar = 0;

        // cache inputs
        this.test = test;
        this.FEssence = FEssence;
        this.FtFinverse = FtFinverse;
        this.perGroupN = perGroupN;
        this.CFixedRand = CFixedRand;
        this.U = U;
        this.thetaNull = thetaNull;
        this.beta = beta;
        this.sigmaError = sigmaError;
        this.sigmaG = sigmaG;

        // calculate intermediate matrices
//        RealMatrix FEssence = params.getDesignEssence().getFullDesignMatrixFixed();
        // TODO: do we ever get here with values that can cause integer overflow,
        //       and if so, does it matter?
        this.N = (double) FEssence.getRowDimension() * perGroupN;
        this.exact = exact;
        try
        {
            // TODO: need to calculate H0, need to adjust H1 for Unirep
            // get design matrix for fixed parameters only
            qF = FEssence.getColumnDimension();
            // a = CFixedRand.getCombinedMatrix().getRowDimension();

            // get fixed contrasts
            RealMatrix Cfixed = CFixedRand.getFixedMatrix();
            RealMatrix CGaussian = CFixedRand.getRandomMatrix();

            // build intermediate terms h1, S
            if (FtFinverse == null)
            {
                FtFinverse = new LUDecomposition(FEssence.transpose().multiply(FEssence)).getSolver().getInverse();
                debug("FEssence", FEssence);
                debug("FtFinverse = (FEssence transpose * FEssence) inverse", FtFinverse);
            } else {
                debug("FtFinverse", FtFinverse);
            }

            RealMatrix PPt = Cfixed.multiply(FtFinverse.scalarMultiply(1/(double) perGroupN)).multiply(Cfixed.transpose());
            debug("Cfixed", Cfixed);
            debug("n = " + perGroupN);
            debug("PPt = Cfixed * FtF inverse * (1/n) * Cfixed transpose", PPt);

            T1 = forceSymmetric(new LUDecomposition(PPt).getSolver().getInverse());
            debug("T1 = PPt inverse", T1);

            FT1 = new CholeskyDecomposition(T1).getL();
            debug("FT1 = Cholesky decomposition (L) of T1", FT1);

            // calculate theta difference
//            RealMatrix thetaNull = params.getTheta();
            RealMatrix C = CFixedRand.getCombinedMatrix();
//            RealMatrix beta = params.getScaledBeta();
//            RealMatrix U = params.getWithinSubjectContrast();

            // thetaHat = C * beta * U
            RealMatrix thetaHat = C.multiply(beta.multiply(U));
            debug("C", C);
            debug("beta", beta);
            debug("U", U);
            debug("thetaHat = C * beta * U", thetaHat);

            // thetaDiff = thetaHat - thetaNull
            RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
            debug("thetaNull", thetaNull);
            debug("thetaDiff = thetaHat - thetaNull", thetaDiff);

            // TODO: specific to HLT or UNIREP
            RealMatrix sigmaStarInverse = getSigmaStarInverse(U, sigmaError, test);
            debug("sigmaStarInverse", sigmaStarInverse);

            RealMatrix H1matrix = thetaDiff.transpose().multiply(T1).multiply(thetaDiff).multiply(sigmaStarInverse);
            debug("H1matrix = thetaDiff transpose * T1 * thetaDiff * sigmaStarInverse", H1matrix);

            H1 = H1matrix.getTrace();
            debug("H1 = " + H1);

            // Matrix which represents the non-centrality parameter as a linear combination of chi-squared r.v.'s.
            S = FT1.transpose().multiply(thetaDiff).multiply(sigmaStarInverse).multiply(thetaDiff.transpose())
                   .multiply(FT1).scalarMultiply(1/H1);
            debug("S = FT1 transpose * thetaDiff * sigmaStar inverse * thetaDiff transpose * FT1 * (1/H1)", S);

            // We use the S matrix to generate the F-critical, numerical df's, and denominator df's
            // for a central F distribution.  The resulting F distribution is used as an approximation
            // for the distribution of the non-centrality parameter.
            // See formulas 18-21 and A8,A10 from Glueck & Muller (2003) for details.
            EigenDecomposition sEigenDecomp = new EigenDecomposition(S);
            sEigenValues = sEigenDecomp.getRealEigenvalues();
            // calculate H0
            if (sEigenValues.length > 0) H0 = H1 * (1 - sEigenValues[0]);
            if (H0 <= 0) H0 = 0;

            // count the # of positive eigen values
            for(double value: sEigenValues)
            {
                if (value > 0) sStar++;
            }
            // TODO: throw error if sStar is <= 0
            // TODO: NO: throw error if sStar != sEigenValues.length instead???
            double stddevG = Math.sqrt(sigmaG.getEntry(0, 0));
            RealMatrix svec = sEigenDecomp.getVT();
            mzSq = svec.multiply(FT1.transpose()).multiply(CGaussian).scalarMultiply(1/stddevG);
            for(int i = 0; i < mzSq.getRowDimension(); i++)
            {
                for (int j = 0; j < mzSq.getColumnDimension(); j++)
                {
                    double entry = mzSq.getEntry(i, j);
                    mzSq.setEntry(i, j, entry*entry); // TODO: is there an apache function to do this?
                }
            }

            debug("exiting initialize normally");
        }
        catch (RuntimeException e)
        {
            LOGGER.warn("exiting initialize abnormally", e);

            throw new PowerException(e.getMessage(),
                    PowerErrorEnum.INVALID_DISTRIBUTION_NONCENTRALITY_PARAMETER);
        }
    }

    /**
     * Reset the total sample size on an existing noncentrality distribution
     * @param perGroupN new per group sample size
     */
    public void setPerGroupSampleSize(int perGroupN)
    throws PowerException {
        initialize(test, FEssence, FtFinverse, perGroupN, CFixedRand, U, thetaNull, beta,
                sigmaError, sigmaG, exact);
    }

    /**
     * Reset the beta matrix on an existing noncentrality distribution
     * @param beta the new beta matrix
     */
    public void setBeta(RealMatrix beta)
    throws PowerException {
        initialize(test, FEssence, FtFinverse, perGroupN, CFixedRand, U, thetaNull, beta,
                sigmaError, sigmaG, exact);
    }

    /**
     * Calculate the probability P(W < w), where W follows the distribution of
     * the non-centrality parameter
     *
     * @param w critical point for which to calculate cumulative probability
     * @return P(W < w)
     */
    public double cdf(double w) throws PowerException
    {
        if (H1 <= 0 || w <= H0) return 0;
        if (H1 - w <= 0) return 1;
        ArrayList<ChiSquareTerm> chiSquareTerms = new ArrayList<ChiSquareTerm>();

        try
        {
            double b0 = 1 - w / H1;
            double m1Positive = 0;
            double m1Negative = 0;
            double m2Positive = 0;
            double m2Negative = 0;
            //
            int numPositive = 0;
            int numNegative = 0;
            double nu;
            double delta;
            double lambda;
            double lastPositiveNoncentrality = 0; // for special cases
            double lastNegativeNoncentrality = 0; // for special cases

            // add in the first chi-squared term in the estimate of the non-centrality
            // (expressed as a sum of weighted chi-squared r.v.s)
            // initial chi-square term is central (delta=0) with N-qf df, and lambda = b0
            nu = N - qF;
            lambda = b0;
            delta = 0;
            chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
            // accumulate terms
            if (lambda > 0)
            {
                // positive terms
                numPositive++;
                lastPositiveNoncentrality = delta;
                m1Positive += lambda * (nu + delta);
                m2Positive += lambda * lambda * 2* (nu + 2*delta);
            }
            else if (lambda < 0)
            {
                // negative terms - we take absolute value of lambda where needed
                numNegative++;
                lastNegativeNoncentrality = delta;
                m1Negative += -1 * lambda * (nu + delta);
                m2Negative += lambda * lambda * 2* (nu + 2*delta);
            }

            // accumulate the remaining terms
            for(int k = 0; k < sStar; k++)
            {
                if (k < sStar)
                {
                    // for k = 1 (well, 0 in java array terms and 1 in the paper) to sStar, chi-square term is
                    // non-central (delta = mz^2), 1 df, lambda = (b0 - kth eigen value of S)
                    nu = 1;
                    lambda = b0 - sEigenValues[k];
                    delta = mzSq.getEntry(k, 0);
                    chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
                }
                else
                {
                    // for k = sStar+1 to a, chi-sqaure term is non-central (delta = mz^2), 1 df,
                    // lambda = b0
                    nu = 1;
                    lambda = b0;
                    delta = mzSq.getEntry(k, 0);
                    chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
                }
                // accumulate terms
                if (lambda > 0)
                {
                    // positive terms
                    numPositive++;
                    lastPositiveNoncentrality = delta;
                    m1Positive += lambda * (nu + delta);
                    m2Positive += lambda * lambda * 2* (nu + 2*delta);
                }
                else if (lambda < 0)
                {
                    // negative terms - we take absolute value of lambda where needed
                    numNegative++;
                    lastNegativeNoncentrality = delta;
                    m1Negative += -1 * lambda * (nu + delta);
                    m2Negative += lambda * lambda * 2* (nu + 2*delta);
                }
                // Note, we deliberately ignore terms for which lambda == 0
            }

            // handle special cases
            if (numNegative == 0) return 0;
            if (numPositive == 0) return 1;

            // special cases
            if (numNegative == 1 && numPositive == 1)
            {
                double Nstar = N - qF + a - 1;
                double Fstar = w / (Nstar * (H1 - w));
                if (lastPositiveNoncentrality >= 0 && lastNegativeNoncentrality == 0)
                {
                    // handle special case: CGaussian = 0, s* = 1
                    NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(Nstar, 1, lastPositiveNoncentrality);
                    return nonCentralFDist.cdf(Fstar);
                }
                else if (lastPositiveNoncentrality == 0 && lastNegativeNoncentrality > 0)
                {
                    // handle special case: CGaussian = 1
                    NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(1,Nstar, lastNegativeNoncentrality);
                    return 1 - nonCentralFDist.cdf(1/Fstar);
                }
            }

            if (exact)
            {
                WeightedSumOfNoncentralChiSquaresDistribution dist    =
                    new WeightedSumOfNoncentralChiSquaresDistribution(chiSquareTerms, ACCURACY);
                return dist.cdf(0);
            }
            else
            {
                // handle general case - Satterthwaite approximation
                double nuStarPositive = 2 * (m1Positive * m1Positive) / m2Positive;
                double nuStarNegative = 2 * (m1Negative * m1Negative) / m2Negative;
                double lambdaStarPositive = m2Positive / (2 * m1Positive);
                double lambdaStarNegative =  m2Negative / (2 * m1Negative);

                // create a central F to approximate the distribution of the non-centrality parameter
                FDistribution centralFDist = new FDistribution(nuStarPositive, nuStarNegative);
                // return power based on the non-central F
                return centralFDist.cumulativeProbability(
                        (nuStarNegative*lambdaStarNegative)/(nuStarPositive*lambdaStarPositive));
            }
        }
        catch (RuntimeException e)
        {
            LOGGER.warn("exiting cdf abnormally", e);

            throw new PowerException(e.getMessage(),
                    PowerErrorEnum.DISTRIBUTION_NONCENTRALITY_PARAMETER_CDF_FAILED);
        }
    }

    /**
     * For this non-centrality distribution, W, this function returns the critical value, w,
     * such that P(W < w).
     *
     * @param probability desired value of P(W < w)
     * @return critical w such that P(W < w)
     */
    public double inverseCDF(double probability)
    {
        if (H1 <= 0) return 0;

        BisectionSolver solver = new BisectionSolver();
        NonCentralityQuantileFunction quantFunc = new NonCentralityQuantileFunction(probability);

        try
        {
            return solver.solve(MAX_ITERATIONS, quantFunc, H0, H1);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to determine non-centrality quantile: " + e.getMessage());
        }
    }

    /**
     * Calculate the inverse of the sigma star matrix
     * @param params GLMM input parameters
     * @return sigma star inverse
     */
    private RealMatrix getSigmaStarInverse(RealMatrix U, RealMatrix sigmaError, Test test)
    {
        // sigma* = U'*sigmaE*U
        RealMatrix sigmaStar = U.transpose().multiply(sigmaError).multiply(U);
        debug("U", U);
        debug("sigmaError", sigmaError);
        debug("sigmaStar = U transpose * sigmaError * U", sigmaStar);

        // TODO: force symmetric?

        if (! MatrixUtils.isPositiveDefinite(sigmaStar)) {
            throw new IllegalArgumentException(NOT_POSITIVE_DEFINITE);
        }

        if (test == Test.HOTELLING_LAWLEY_TRACE)
        {
            return new LUDecomposition(sigmaStar).getSolver().getInverse();
        }
        else
        {
            // stat should only be UNIREP (uncorrected, box, GG, or HF) at this point
            // (exception is thrown by valdiateParams otherwise)
            int b = sigmaStar.getColumnDimension();
            // get discrepancy from sphericity for unirep test
            double sigmaStarTrace = sigmaStar.getTrace();
            double sigmaStarSquaredTrace = sigmaStar.multiply(sigmaStar).getTrace();
            double epsilon = (sigmaStarTrace*sigmaStarTrace) / ((double) b * sigmaStarSquaredTrace);
            RealMatrix identity = org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(b);
            return identity.scalarMultiply((double) b * epsilon / sigmaStarTrace);
        }
    }

    /**
     * Get the upper integral bound
     * @return H1
     */
    public double getH1()
    {
        return H1;
    }

    /**
     * Get the lower integral bound
     * @return H0
     */
    public double getH0()
    {
        return H0;
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
     * A convenience method for DEBUG logging of a message
     * and a throwable.
     *
     * @param message The message.
     * @param t       The throwable.
     */
    private static void debug(Object message, Throwable t) {
        LOGGER.debug(message, t);
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
