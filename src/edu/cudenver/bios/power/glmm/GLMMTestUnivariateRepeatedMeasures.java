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
import java.util.Arrays;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.Precision;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.GramSchmidtOrthonormalization;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;

/**
 * Implementation of the uncorreected univariate approach to repeated measures test
 * (UNIREP) for the general linear multivariate model.  This is also the base class
 * for the corrected tests (Box, Geisser-Greenhouse, Huynh-Feldt)
 *
 * @author Sarah Kreidler
 *
 */
public class GLMMTestUnivariateRepeatedMeasures extends GLMMTest
{
    // if sigma is estimated, then this value will be set to totalN_est - rank_est
    // where totalN_est is the sample size for the data set from which sigma was estimated
    // and rank_est is the rank of the design matrix in the data set used for estimation
    protected int nuEst =  0;
    // cache some of the eigenvalue information for use in the GG and Huynh-Feldt
    protected ArrayList<EigenValueMultiplicityPair> distinctSigmaStarEigenValues = new ArrayList<EigenValueMultiplicityPair>();
    protected double[] sigmaStarEigenValues = null;
    protected double sumLambda = 0;
    protected double sumLambdaSquared = 0;
    protected int rankC = 0; // rank of C
    protected int rankU = 0; // rank of U

    private double eigenTolerance = 1.0E-15;

    /* correction factors for sphericity
     * these differ for data analysis, power under the null and power under the alternative
     * also, the values depend on whether sigma is known or estimated
     * For full details please see
     *
     * Gribbin MJ (2007). Better Power Methods for the Univariate Approach to Repeated Measures.
     * Ph.D. thesis, University of North Carolina at Chapel Hill.
     * (see p61, equation 50-52 and table 4.1)
     *
     * Per Gribbin, epsilonD is the epsilon value under the null case, and epsilonN
     * is the sphericity parameter in the non-null case
     */

    //  CDF approximation method;
    protected UnivariateCdfApproximation cdfMethod;
    // epsilon approximation information - this only applies to HF, GG, but
    // added to this class to avoid duplicated code
    protected UnivariateEpsilonApproximation epsilonMethod;

    // sphericity components when sigma known
    protected double epsilonD = Double.NaN;
    protected double epsilonN = Double.NaN;
    // sphericity components when sigma estimated
    protected double epsilonTildeN = Double.NaN;
    protected double epsilonTildeR = Double.NaN;
    protected double dataAnalysisNDFCorrection = 1;
    protected double dataAnalysisDDFCorrection = 1;
    protected double powerNullNDFCorrection = 1;
    protected double powerNullDDFCorrection = 1;
    protected double powerAlternativeNDFCorrection = 1;
    protected double powerAlternativeDDFCorrection = 1;
    protected double noncentralityCorrection = 1;

    // class for tracking eigen value information during epsilon calculation
    protected class EigenValueMultiplicityPair
    {
        public double eigenValue;
        public double multiplicity;

        public EigenValueMultiplicityPair(double eigenValue, double multiplicity)
        {
            this.eigenValue = eigenValue;
            this.multiplicity = multiplicity;
        }
    };

    /**
     * Create a UNIREP test object for the specified parameters
     * @param params GLMM input parameters
     */
    public GLMMTestUnivariateRepeatedMeasures(FApproximation fMethod,
            UnivariateCdfApproximation cdfMethod, UnivariateEpsilonApproximation epsilonMethod,
            RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
            FixedRandomMatrix C, RealMatrix U, RealMatrix thetaNull,
            RealMatrix beta, RealMatrix sigmaError, int nuEst)
    {
        super(fMethod, Xessence, XtXInverse, perGroupN, rank,
                C, U, thetaNull, beta, sigmaError);
        this.cdfMethod = cdfMethod;
        this.epsilonMethod = epsilonMethod;
        this.nuEst = nuEst;
        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();

        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateEpsilon();
        // calculate the adjustment factors for degrees of freedom, noncentrality
        // as described in Gribbin.  Note that calculateEpsilon must be called before
        // these functions
        calculateNDFCorrection();
        calculateDDFCorrection();
        calculateNoncentralityCorrection();
    }

    /**
     * Create a UNIREP test object for data analysis
     * @param params GLMM input parameters
     */
    public GLMMTestUnivariateRepeatedMeasures(FApproximation fMethod,
            UnivariateCdfApproximation cdfMethod, UnivariateEpsilonApproximation epsilonMethod,
            RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
            RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        super(fMethod, X, XtXInverse, rank, Y,  C, U, thetaNull);
        this.cdfMethod = cdfMethod;
        this.epsilonMethod = epsilonMethod;
        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();

        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateEpsilon();
        // calculate the adjustment factors for degrees of freedom, noncentrality
        // as described in Gribbin.  Note that calculateEpsilon must be called before
        // these functions
        calculateNDFCorrection();
        calculateDDFCorrection();
        calculateNoncentralityCorrection();
    }

    @Override
    public void setPerGroupSampleSize(int perGroupN)
    {
        super.setPerGroupSampleSize(perGroupN);
        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateEpsilon();
        // calculate the adjustment factors for degrees of freedom, noncentrality
        // as described in Gribbin.  Note that calculateEpsilon must be called before
        // these functions
        calculateNDFCorrection();
        calculateDDFCorrection();
        calculateNoncentralityCorrection();
    }

    /**
     * Calculate the denominator degrees of freedom for the UNIREP, based on
     * whether the null or alternative hypothesis is assumed true.
     *
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();

        double df = b*(totalN - rank);
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected test, we adjust by epsilon only for
        // power analysis under the alternative.  The ddf are the same for power
        // under the null and for data analysis
        switch (type)
        {
        case POWER_ALTERNATIVE:
            df *= powerAlternativeDDFCorrection;
            break;
        case POWER_NULL:
            df *= powerNullDDFCorrection;
            break;
        case DATA_ANALYSIS_NULL:
            df *= dataAnalysisDDFCorrection;
            break;
        }

        return df;
    }

    /**
     * Calculate the non-centrality parameter for the UNIREP, based on
     * whether the null or alternative hypothesis is assumed true.
     *
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    @Override
    public double getNonCentrality(DistributionType type)
    {
        //double a = C.getRowDimension();
        double b = U.getColumnDimension();
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares();
        // TODO: cache sig star, lam bar
        RealMatrix sigmaStar = U.transpose().multiply(sigmaError.multiply(U));
        double lambdaBar = sigmaStar.getTrace() / rankU;

        return (hypothesisSumOfSquares.getTrace() / (lambdaBar / noncentralityCorrection));
        // calculate non-centrality and adjust for sphericity

        //return a*b*getObservedF(type)*(nuEst > 0 ? this.estimatedSigmaCorrection : this.unirepEpsilon);
    }

    /**
     * Calculate the numerator degrees of freedom for the UNIREP, based on
     * whether the null or alternative hypothesis is assumed true.
     *
     * @param type distribution type
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = C.getRowDimension();
        double b = U.getColumnDimension();

        double df = a * b;
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected, we adjust by epsilon only for
        // power analysis under the alternative.  The ndf are the same for power
        // under the null and for data analysis
        switch (type)
        {
        case POWER_ALTERNATIVE:
            df *= powerAlternativeNDFCorrection;
            break;
        case POWER_NULL:
            df *= powerNullNDFCorrection;
            break;
        case DATA_ANALYSIS_NULL:
            df *= dataAnalysisNDFCorrection;
            break;
        }

        return df;
    }

    /**
     * Calculate the observed F for the UNIREP, based on
     * whether the null or alternative hypothesis is assumed true.
     *
     * @param type distribution type
     * @return observed F
     * @throws IllegalArgumentException
     */
    @Override
    public double getObservedF(DistributionType type)
    {
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares();
        RealMatrix errorSumOfSquares = getErrorSumOfSquares();

        double association = 0.0;

        double REP = getUnirep(hypothesisSumOfSquares, errorSumOfSquares);
        association = (REP / (1 + REP));

        double ddf = getDenominatorDF(type);
        double ndf = getNumeratorDF(type);
        return ((association) / ndf) / ((1 - association) / ddf);
    }

    /**
     * Compute a Univariate Approach to Repeated Measures statistic
     *
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    protected double getUnirep(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new IllegalArgumentException("Failed to compute Unirep statistic: hypothesis and error matrices must be square and same dimensions");

        return H.getTrace() / E.getTrace();
    }

    /**
     * Calculate the epsilon to correct for violations of sphericity
     */
    protected void calculateEpsilon()
    {
        distinctSigmaStarEigenValues.clear();
        rankC = new SingularValueDecomposition(C).getRank();
        rankU = new SingularValueDecomposition(U).getRank();

        // get the sigmaStar matrix: U' *sigmaError * U
       RealMatrix  sigmaStar = U.transpose().multiply(sigmaError.multiply(U));

        // ensure symmetry
        sigmaStar = sigmaStar.add(sigmaStar.transpose()).scalarMultiply(0.5);

        // get the eigen values of the normalized sigmaStar matrix
        sigmaStarEigenValues =
            new EigenDecomposition(sigmaStar.scalarMultiply(1/sigmaStar.getTrace())).getRealEigenvalues();
        if (sigmaStarEigenValues.length <= 0) throw new IllegalArgumentException("Failed to compute eigenvalues for sigma* matrix");
        Arrays.sort(sigmaStarEigenValues);
        // get the trace of sigma* and sigma* squared
        // to avoid looping over the eigenvalues twice, we also calculate the multiplicity for distinct eigenvalues

        // initialize values for the first eigen value
        double first = sigmaStarEigenValues[0];
        distinctSigmaStarEigenValues.add(new EigenValueMultiplicityPair(first, 1));
        sumLambda = first; // sum of the eigenvalues (i.e. trace of sigmaStar)
        sumLambdaSquared = first * first;  // sum of squared eignevalues (trace of sigmaStar squared)

        // loop over remaining eigen values, saving distinct eigen values
        for(int i = 1; i < sigmaStarEigenValues.length; i++)
        {
            double value = sigmaStarEigenValues[i];
            // build the sum & sum of squares of eigen values
            sumLambda += value;
            sumLambdaSquared += value * value;

            // determine if this is a distinct eigen value and calculate multiplicity
            EigenValueMultiplicityPair prev = distinctSigmaStarEigenValues.get(distinctSigmaStarEigenValues.size()-1);
            if (Math.abs(prev.eigenValue - value) > eigenTolerance)
            {
                // found new distinct eigen value
                distinctSigmaStarEigenValues.add(new EigenValueMultiplicityPair(value, 1));
            }
            else
            {
                // repeat of same eigenvalue, so  increment the multiplicity
                prev.multiplicity++;
            }
        }

        //get the hypothesis sum of squares times 1/a
        RealMatrix H  = getHypothesisSumOfSquares();
        RealMatrix HDivA = H.scalarMultiply((double)1/(double)rankC);

        double sigStarTrace = sigmaStar.getTrace();
        double sigStarSqTrace = MatrixUtils.getSumOfSquares(sigmaStar); //sigmaStar.transpose().multiply(sigmaStar).getTrace();
        RealMatrix sigStarH = sigmaStar.multiply(getHypothesisSumOfSquares());
        double sigStarHTrace = sigStarH.getTrace();

        // calculate estimate of epsilon (correction for violation of spehericity assumption)
        epsilonD = (sigStarTrace*sigStarTrace) / (rankU * (sigStarSqTrace));
        epsilonN = (sigStarTrace * sigStarTrace +
                2 * sigStarTrace * HDivA.getTrace()) /
                (rankU * (sigStarSqTrace +
                        2 * (sigmaStar.multiply(HDivA).getTrace())));
        if (cdfMethod == UnivariateCdfApproximation.MULLER_BARTON_APPROX)
            epsilonN = epsilonD;

        if (nuEst > 0)
        {

//               EPSNHAT_NUM =  NU_EST # (Q3 # NU_EST # (NU_EST+{1}) + Q1 # Q2 # ({2}#MULT/A) - Q4 # {2} # NU_EST) ;
//               EPSNHAT_DEN = (Q4 # NU_EST # NU_EST + Q5 #(2#MULT/A) - Q3# NU_EST)# NU_EST ;
//               EPSNHAT = EPSNHAT_NUM / (B # EPSNHAT_DEN);
            double multiplier = nuEst * nuEst + nuEst - 2;
            double numerator = nuEst * (nuEst*(nuEst + 1)*sigStarTrace*sigStarTrace -
                    2 * nuEst * sigStarSqTrace +
                    2 * multiplier * H.getTrace() * sigStarTrace / rankC);
            double denominator = nuEst * (nuEst*nuEst*sigStarSqTrace - nuEst *sigStarTrace*sigStarTrace +
                    2 *multiplier * sigStarHTrace / rankC);

            epsilonTildeN = numerator / ( rankU *denominator);

            epsilonTildeR = ((nuEst + 1)*rankU*epsilonD - 2) / (rankU*(nuEst - rankU * epsilonD));
            epsilonTildeR = (epsilonTildeR < 1 ? epsilonTildeR : 1);
        }
    }

    /**
     * Ensure that the within subject contrast is orthonormal for all
     * UNIREP tests
     */
    protected void createOrthonormalU()
    {
        RealMatrix UtU = U.transpose().multiply(U);
        double upperLeft = UtU.getEntry(0, 0);
        if (upperLeft != 0) UtU = UtU.scalarMultiply(1/upperLeft);

        RealMatrix diffFromIdentity =
            UtU.subtract(org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(UtU.getRowDimension()));

        // get maximum value in U'U
        double maxValue = Double.NEGATIVE_INFINITY;
        for(int r = 0; r < diffFromIdentity.getRowDimension(); r++)
        {
            for(int c = 0; c < diffFromIdentity.getColumnDimension(); c++)
            {
                double entryVal = diffFromIdentity.getEntry(r, c);
                if (entryVal > maxValue) maxValue = entryVal;
            }
        }

        if (maxValue > Precision.SAFE_MIN)
        {
            // U matrix deviates from Identity, so create one that is orthonormal
            RealMatrix Utmp = new GramSchmidtOrthonormalization(U).getQ();
            U = Utmp;
        }
    }

    /**
     * Calculate the correction factors for numerator degrees of
     * freedom for data analysis, power under the null and power
     * under the alternative
     */
    protected void calculateNDFCorrection()
    {
        dataAnalysisNDFCorrection = 1;
        powerNullNDFCorrection = 1;
        if (nuEst <= 0)
            powerAlternativeNDFCorrection = epsilonN;
        else
            powerAlternativeNDFCorrection = epsilonTildeN;
    }

    /**
     * Calculate the correction factors for denominator degrees of
     * freedom for data analysis, power under the null and power
     * under the alternative
     */
    protected void calculateDDFCorrection()
    {
        dataAnalysisDDFCorrection = 1;
        powerNullDDFCorrection = 1;
        powerAlternativeDDFCorrection = epsilonD;
    }

    /**
     * Calculate the correction factors for noncentrality
     * parameter.  This is only relevant for power under the alternative.
     */
    protected void calculateNoncentralityCorrection()
    {
        if (nuEst <= 0)
            noncentralityCorrection = epsilonN;
        else
            noncentralityCorrection = epsilonTildeN;
    }

    /**
     * Get the degrees of freedom for the Chi-sqaured distribution used
     * in building confidence limits on power as described by Gribbin.
     *
     * @return
     */
    public double getConfidenceLimitsDegreesOfFreedom()
    {
        return rankU * nuEst * powerAlternativeDDFCorrection / noncentralityCorrection;
    }
}
