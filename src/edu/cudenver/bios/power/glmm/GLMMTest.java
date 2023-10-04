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

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;

/**
 * Abstract base class for statistical tests for the GLMM
 * @author Sarah Kreidler
 *
 */
public abstract class GLMMTest
{
    private static final String OVERFLOW_ERROR_MESSAGE =
        "Sorry, an overflow occurred. Please try reducing the smallest group size.";

    // type of approximation to use for unirep
    public enum UnivariateCdfApproximation
    {
        MULLER_BARTON_APPROX,
        MULLER_EDWARDS_TAYLOR_APPROX,
        MULLER_EDWARDS_TAYLOR_EXACT,
        MULLER_EDWARDS_TAYLOR_EXACT_APPROX
    };

    // type of approximation to use for mean epsilon
    // for Huynh-Feldt and Geisser-Greenhouse
    public enum UnivariateEpsilonApproximation
    {
        MULLER_BARTON_APPROX,
        MULLER_EDWARDS_TAYLOR_APPROX,
    };

    // available f approximation methods
    public enum FApproximation
    {
        NONE,
        PILLAI_ONE_MOMENT,
        PILLAI_ONE_MOMENT_OMEGA_MULT,
        MCKEON_TWO_MOMENT,
        MCKEON_TWO_MOMENT_OMEGA_MULT,
        MULLER_TWO_MOMENT,
        MULLER_TWO_MOMENT_OMEGA_MULT,
        RAO_TWO_MOMENT,
        RAO_TWO_MOMENT_OMEGA_MULT
    };

    // for the unirep test, the degrees of freedom change depending on two
    // factors - data analysis (i.e. simulation or model fit) vs. power analysis,
    // and whether we need df for the F distribution under the null or the
    // alternative hypothesis
    public enum DistributionType
    {
        DATA_ANALYSIS_NULL,
        POWER_NULL,
        POWER_ALTERNATIVE
    };

    /**
     * container class for model fit
     */
    public class ModelFit
    {
        public RealMatrix sigma;
        public RealMatrix beta;
        public double numeratorDF;
        public double denominatorDF;
        public double Fvalue;
        public double Pvalue;

        public ModelFit(double Pvalue, double Fvalue,
                double numeratorDF, double denominatorDF,
                RealMatrix sigma, RealMatrix beta)
        {
            this.sigma = sigma;
            this.beta = beta;
            this.Pvalue = Pvalue;
            this.Fvalue = Fvalue;
            this.numeratorDF = numeratorDF;
            this.denominatorDF = denominatorDF;
        }
    }

    // approximation information
    protected FApproximation fMethod;

    // study design matrices
    protected RealMatrix Xessence;
    protected RealMatrix XtXInverse;
    protected double totalN; // total sample size
    protected double rank; // rank of the design matrix
    // contrasts
    // we hang onto to the fixed portion of the C matrix to allow updating the
    // per group N when we have a GLMM(F,g) design
    protected RealMatrix CFixed;
    protected RealMatrix C; // between subject
    protected RealMatrix U; // within subject
    // null hypothesis values
    protected RealMatrix thetaNull;
    // estimated beta parameters - only for power analysis
    protected RealMatrix beta;
    // estimated sigma error matrix
    protected RealMatrix sigmaError;
    // the M matrix: [C(X'X)-1C']-1
    protected RealMatrix M = null;
    // flag indicating whether the design is multivariate or univariate
    protected boolean multivariate = false;
    /**
     * Create a statistical test for to compute power
     * @param params GLMM input parameters
     */
    public GLMMTest(FApproximation fMethod,
            RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
            FixedRandomMatrix C, RealMatrix U, RealMatrix thetaNull,
            RealMatrix beta, RealMatrix sigmaError)
    {
        this.fMethod = fMethod;
        this.Xessence = Xessence;
        this.XtXInverse = XtXInverse;
        int xrd = Xessence.getRowDimension();
        if (xrd != 0 && perGroupN > Integer.MAX_VALUE/xrd) {
            throw new GLMMTestException(OVERFLOW_ERROR_MESSAGE);
        }
        this.totalN = (double) xrd * perGroupN;
        this.rank = rank;
        this.CFixed = C.getFixedMatrix();
        this.C = C.getCombinedMatrix();
        this.U = U;
        this.thetaNull = thetaNull;
        this.beta = beta;
        this.sigmaError = sigmaError;

        // check if uni/multivariate design
        multivariate = (beta.getColumnDimension() > 1);

        // cache the value of M
        RealMatrix CFixed = C.getFixedMatrix();
        RealMatrix cxxcEssence = CFixed.multiply((XtXInverse).multiply(CFixed.transpose()));
        RealMatrix cxxcEssenceInverse = new LUDecomposition(cxxcEssence).getSolver().getInverse();
        this.M = cxxcEssenceInverse.scalarMultiply(perGroupN);
    }

    /**
     * Create a test to perform data analysis
     * @param fMethod
     * @param cdfMethod
     * @param X
     * @param XtXInverse
     * @param rank
     * @param Y
     * @param C
     * @param U
     * @param thetaNull
     */
    public  GLMMTest(FApproximation fMethod,
            RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
            RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        this.fMethod = fMethod;
        this.Xessence = null;
        this.XtXInverse = XtXInverse;
        if (this.XtXInverse == null)
        {
            this.XtXInverse =
                new LUDecomposition(X.transpose().multiply(X)).getSolver().getInverse();
        }
        this.totalN = X.getRowDimension();
        this.rank = rank;
        this.C = C;
        this.U = U;
        this.thetaNull = thetaNull;
        this.beta = this.XtXInverse.multiply(X.transpose()).multiply(Y);
        RealMatrix YHat = X.multiply(this.beta);
        RealMatrix Ydiff = Y.subtract(YHat);
        this.sigmaError = (Ydiff.transpose().multiply(Ydiff)).scalarMultiply(((double) 1/(double) (totalN - rank)));

        // check if uni/multivariate design
        multivariate = (Y.getColumnDimension() > 1);

        // cache the value of M
        RealMatrix cxxcEssence = C.multiply((this.XtXInverse).multiply(C.transpose()));
        M = new LUDecomposition(cxxcEssence).getSolver().getInverse();
    }

    /**
     * Reset the total sample size for this test
     * @param totalN total sample size
     */
    public void setPerGroupSampleSize(int perGroupN)
    {
        int xrd = Xessence.getRowDimension();
        // TODO: do we ever get here with values that can cause integer overflow,
        //       and if so, does it matter?
        /*
        if (xrd != 0 && perGroupN > Integer.MAX_VALUE/xrd) {
            throw new GLMMTestException(OVERFLOW_ERROR_MESSAGE);
        }
        */
        this.totalN = (double) xrd * perGroupN;
        RealMatrix cxxcEssence = null;
        if (CFixed != null)
            cxxcEssence = CFixed.multiply((XtXInverse).multiply(CFixed.transpose()));
        else
            cxxcEssence = C.multiply((XtXInverse).multiply(C.transpose()));
        RealMatrix cxxcEssenceInverse = new LUDecomposition(cxxcEssence).getSolver().getInverse();
        this.M = cxxcEssenceInverse.scalarMultiply(perGroupN);
    }

    /**
     * Reset the total sample size for this test
     * @param totalN total sample size
     */
    public void setBeta(RealMatrix beta)
    {
        this.beta = beta;
    }

    /**
     * Fit the model
     * @return model fit object
     */
    public ModelFit getModelFit()
    {
        double fobs = getObservedF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

        // get the p-value from a central F distribution
        double ndf = getNumeratorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
        if (Double.isNaN(ndf)) {
            throw new IllegalArgumentException("numerator DF is NaN");
        }
        double ddf = getDenominatorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
        if (Double.isNaN(ddf)) {
            throw new IllegalArgumentException("denominator DF is NaN");
        }

        FDistribution fdist = new FDistribution(ndf, ddf);
        double pvalue = 1 - fdist.cumulativeProbability(fobs);

           return new ModelFit(pvalue, fobs, ndf, ddf,
                sigmaError, beta);
    }

    /**
     * Calculate the critical F value under the specified distribution
     *
     * @param type distribution type
     * @param alpha type I error level
     * @return critical F
     *
     */
    public double getCriticalF(DistributionType type, double alpha)
    throws IllegalArgumentException
    {
        double ndf = getNumeratorDF(type);
        if (Double.isNaN(ndf)) {
            throw new IllegalArgumentException("numerator DF is NaN");
        }
        double ddf = getDenominatorDF(type);
        if (Double.isNaN(ddf)) {
            throw new IllegalArgumentException("denominator DF is NaN");
        }

        FDistribution centralFDist = new FDistribution(ndf, ddf);
        double fcrit = centralFDist.inverseCumulativeProbability(1 - alpha);
        return fcrit;
    }

    /**
     * Prototype for getting numerator degrees of freedom for the
     * statistical test under the specified distribution type.
     *
     * @param type distribution type
     * @return numerator degrees of freedom
     */
    abstract public double getNumeratorDF(DistributionType type);

    /**
     * Prototype for getting denominator degrees of freedom for the
     * statistical test under the specified distribution type.
     *
     * @param type distribution type
     * @return denominator degrees of freedom
     */
    abstract public double getDenominatorDF(DistributionType type);

    /**
     * Prototype for getting observed F value for the
     * statistical test under the specified distribution type.
     *
     * @param type distribution type
     * @return observed F
     */
    abstract public double getObservedF(DistributionType type);

    /**
     * Prototype for getting the non-centrality parameter for the
     * statistical test under the specified distribution type.
     *
     * @param type distribution type
     * @return non-centrality parameter
     */
    abstract public double getNonCentrality(DistributionType type);

    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    public RealMatrix getHypothesisSumOfSquares()
    {
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(beta.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(thetaNull);

        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(M.multiply(thetaDiff));

        return hss;

    }

    /**
     * Calculate the sum of squares error matrix (the E matrix)
     *
     * @param params matrices input by the user
     * @return error sum of squares
     */
    public RealMatrix getErrorSumOfSquares()
    {
        return U.transpose().multiply(sigmaError.multiply(U)).scalarMultiply(totalN-rank);
    }

    /**
     * Indicates if the test is a univariate or multivariate design
     * based on the columns of beta.  Only available for
     * power analysis
     */
    public boolean isMultivariate()
    {
        return multivariate;
    }

    /**
     * Get the beta matrix
     * @return beta matrix
     */
    public RealMatrix getBeta() {
        return beta;
    }

    /**
     * A runtime exception this class may throw.
     */
    public static final class GLMMTestException extends RuntimeException {
        /**
         * Construct an instance of this class.
         *
         * @param message The exception message.
         */
        public GLMMTestException(String message) {
            super(message);
        }
    }
}
