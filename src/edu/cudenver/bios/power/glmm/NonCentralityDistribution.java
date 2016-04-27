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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.distribution.ChiSquareTerm;
import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.distribution.WeightedSumOfNoncentralChiSquaresDistribution;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.PowerErrorEnum;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;

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
    private static final int MAX_ITERATIONS = Integer.MAX_VALUE;
    private static final double TOLERANCE = 0.000000000001;
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
    int N;
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
                return Double.NaN;
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
        this.N = FEssence.getRowDimension() * perGroupN;
        this.exact = exact;
        try
        {
            // TODO: need to calculate H0, need to adjust H1 for Unirep
            // get design matrix for fixed parameters only
            qF = FEssence.getColumnDimension();
            a = CFixedRand.getCombinedMatrix().getRowDimension();
            // get fixed contrasts
            RealMatrix Cfixed = CFixedRand.getFixedMatrix();
            RealMatrix CGaussian = CFixedRand.getRandomMatrix();
            // build intermediate terms h1, S
            if (FtFinverse == null)
            {
                FtFinverse = new LUDecomposition(FEssence.transpose().multiply(FEssence)).getSolver().getInverse();
            }
            //CF*FPFINV*CF`
            RealMatrix PPt = Cfixed.multiply(FtFinverse.scalarMultiply(1/(double) perGroupN)).multiply(Cfixed.transpose());
            T1 = new LUDecomposition(PPt).getSolver().getInverse();
            FT1 = new CholeskyDecomposition(T1).getL();
            // calculate theta difference
//            RealMatrix thetaNull = params.getTheta();
            RealMatrix C = CFixedRand.getCombinedMatrix();
//            RealMatrix beta = params.getScaledBeta();
//            RealMatrix U = params.getWithinSubjectContrast();
            // thetaHat = C * beta * U
            RealMatrix thetaHat = C.multiply(beta.multiply(U));
            // thetaHat - thetaNull.
            RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
            // TODO: specific to HLT or UNIREP
            RealMatrix sigmaStarInverse = getSigmaStarInverse(U, sigmaError, test);
            RealMatrix H1matrix = thetaDiff.transpose().multiply(T1).multiply(thetaDiff).multiply(sigmaStarInverse);
            H1 = H1matrix.getTrace();
            // matrix which represents the non-centrality parameter as a linear combination of chi-squared r.v.'s
            S = FT1.transpose().multiply(thetaDiff).multiply(sigmaStarInverse).multiply(thetaDiff.transpose()).multiply(FT1).scalarMultiply(1/H1);
            // we use the S matrix to generate the F-critical, numerical df's, and denominator df's
            // for a central F distribution.  The resulting F distribution is used as an approximation
            // for the distribution of the non-centrality parameter
            // See formulas 18-21 and A8,A10 from Glueck & Muller (2003) for details
            EigenDecomposition sEigenDecomp = new EigenDecomposition(S, TOLERANCE);
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
        }
        catch (Exception e)
        {
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
        catch (Exception e)
        {
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
            RealMatrix identity = MatrixUtils.createRealIdentityMatrix(b);
            return identity.scalarMultiply((double) b * epsilon / sigmaStarTrace);
        }
    }

    /**
     * Get the upper intergral bound
     * @return H1
     */
    public double getH1()
    {
        return H1;
    }

    /**
     * Get the lower intergral bound
     * @return H0
     */
    public double getH0()
    {
        return H0;
    }
}
