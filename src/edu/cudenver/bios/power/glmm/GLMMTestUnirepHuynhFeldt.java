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

import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;

/**
 * Implementation of the univariate approach to repeated measures test
 * with Huynh-Feldt correction (UNIREP-HF) for the general linear multivariate model.
 *
 * @see GLMMTestUnivariateRepeatedMeasures
 * @author Sarah Kreidler
 *
 */
public class GLMMTestUnirepHuynhFeldt extends GLMMTestUnivariateRepeatedMeasures
{
    protected static final double TOLERANCE = 0.000000000001;
    double expectedEpsilon = Double.NaN;

    /**
     * Create a UNIREP-HF test object for the specified parameters
     * @param params GLMM input parameters
     */
    public GLMMTestUnirepHuynhFeldt(FApproximation fMethod,
            UnivariateCdfApproximation cdfMethod, UnivariateEpsilonApproximation epsilonMethod,
            RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
            FixedRandomMatrix C, RealMatrix U, RealMatrix thetaNull,
            RealMatrix beta, RealMatrix sigmaError, int nuEst)
    {
        super(fMethod, cdfMethod, epsilonMethod, Xessence, XtXInverse, perGroupN, rank,
                C, U, thetaNull, beta, sigmaError, nuEst);
    }

    /**
     * Create a UNIREP test object for data analysis
     * @param params GLMM input parameters
     */
    public GLMMTestUnirepHuynhFeldt(FApproximation fMethod,
            UnivariateCdfApproximation cdfMethod, UnivariateEpsilonApproximation epsilonMethod,
            RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
            RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        super(fMethod, cdfMethod, epsilonMethod, X, XtXInverse, rank, Y,  C, U, thetaNull);
    }

    /**
     * Calculate the Huynh-Feldt epsilon to correct for violations
     * of sphericity
     */
    @Override
    protected void calculateEpsilon()
    {
        super.calculateEpsilon();
        double b = rankU;

        if (this.epsilonMethod == UnivariateEpsilonApproximation.MULLER_BARTON_APPROX)
        {
            // calculate the expected value of the epsilon estimate
            // E[h(lambda)] = h(lambda) + g1 / (N - r)
            // h(lambda) = h1(lambda) / (b*h2(lambda)
            // see Muller, Barton (1989) for details
            double h1 = totalN * sumLambda * sumLambda - 2 * sumLambdaSquared;
            double h2 = (totalN - rank) * sumLambdaSquared - (sumLambda * sumLambda);
            double g1 = 0;
            for(int i = 0; i < distinctSigmaStarEigenValues.size(); i++)
            {
                EigenValueMultiplicityPair evmI = distinctSigmaStarEigenValues.get(i);
                // derivatives of sub-equations comprising epsilon estimator
                double h1firstDerivative = (2 * totalN * sumLambda) - (4 *evmI.eigenValue);
                double h1secondDerivative = 2 * totalN - 4;

                double h2firstDerivative = (2 * (totalN - rank) * evmI.eigenValue) - (2 * sumLambda);
                double h2secondDerivative = 2 * (totalN - rank) - 2;

                // derivatives of estimate of epsilon
                double firstDerivative = ((h1firstDerivative) - ((h1 * h2firstDerivative) / h2)) / (h2 *b);
                double secondDerivative =
                    (h1secondDerivative - ((2*h1firstDerivative*h2firstDerivative)/(h2)) +
                            (2 * h1 * h2firstDerivative * h2firstDerivative)/(h2*h2) -
                            (h1*h2secondDerivative)/(h2)) /
                            (h2*b);

                // accumulate the first term of g1 (sum over distinct eigen vals of 1st derivative * eigen val ^2 * multiplicity)
                g1 += secondDerivative * evmI.eigenValue * evmI.eigenValue * evmI.multiplicity;
                // loop over elements not equal to current
                for(int j = 0; j < distinctSigmaStarEigenValues.size(); j++)
                {
                    if (i != j)
                    {
                        EigenValueMultiplicityPair evmJ = distinctSigmaStarEigenValues.get(j);
                        // accumulate second term of g1
                        g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
                                (evmI.eigenValue - evmJ.eigenValue));
                    }
                }
            }

            expectedEpsilon = (totalN*b*epsilonD - 2)/(b*(totalN - rank - b*epsilonD))  + g1 / (totalN - rank);
        }
        else
        {
            double sum = 0;
            for(int i = 0; i < sigmaStarEigenValues.length; i++)
            {
                for(int j = 0; j < sigmaStarEigenValues.length; j++)
                {
                    sum += sigmaStarEigenValues[i] * sigmaStarEigenValues[j];
                }
            }

            double nu = totalN - rank;
            double expT1 = (2 * nu * sumLambdaSquared) + (nu * nu * sumLambda * sumLambda);
            double expT2 = ((nu * (nu + 1)) * sumLambdaSquared) + (nu * sum);

            double numerator = (1/b) * (((nu + 1)*expT1) - 2 * expT2);
            double denominator = (nu * expT2) - expT1;

            expectedEpsilon = numerator / denominator;
        }

        // ensure that expected value is within bounds 1/b to 1
        if (!Double.isNaN(expectedEpsilon))
        {
            if (expectedEpsilon < 1/b)
            {
                expectedEpsilon = 1/b;
            }
            else if (expectedEpsilon > 1)
            {
                expectedEpsilon = 1;
            }
        }

    }

    /**
     * Calculate the correction factors for numerator degrees of
     * freedom for data analysis, power under the null and power
     * under the alternative
     */
    @Override
    protected void calculateNDFCorrection()
    {
        dataAnalysisNDFCorrection = epsilonD;
        if (nuEst <= 0)
        {
            powerNullNDFCorrection = expectedEpsilon;
            powerAlternativeNDFCorrection = epsilonN;
        }
        else
        {
            // sigma estimated
            powerNullNDFCorrection = epsilonTildeR;
            powerAlternativeNDFCorrection = epsilonTildeN;
        }
    }

    /**
     * Calculate the correction factors for denominator degrees of
     * freedom for data analysis, power under the null and power
     * under the alternative
     */
    @Override
    protected void calculateDDFCorrection()
    {
        dataAnalysisDDFCorrection = epsilonD;
           powerAlternativeDDFCorrection = epsilonD;
        if (nuEst <= 0)
        {
            powerNullDDFCorrection = expectedEpsilon;
        }
        else
        {
            // sigma estimated
            powerNullDDFCorrection = epsilonTildeR;
        }
    }

    /**
     * Calculate the correction factors for noncentrality
     * parameter.  This is only relevant for power under the alternative.
     */
    @Override
    protected void calculateNoncentralityCorrection()
    {
        if (nuEst <= 0)
            noncentralityCorrection = epsilonN;
        else
            noncentralityCorrection = epsilonTildeN;
    }

}
