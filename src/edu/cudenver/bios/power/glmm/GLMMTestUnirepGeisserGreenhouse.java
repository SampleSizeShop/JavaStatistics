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

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;

import edu.cudenver.bios.matrix.FixedRandomMatrix;

/**
 * Implementation of the univariate approach to repeated measures test
 * with Geisser-Greenhouse correction (UNIREP-GG) for the general linear multivariate model.
 *
 * @see GLMMTestUnivariateRepeatedMeasures
 * @author Sarah Kreidler
 *
 */
public class GLMMTestUnirepGeisserGreenhouse extends GLMMTestUnivariateRepeatedMeasures
{
    protected static final double TOLERANCE = 0.000000000001;
    private double expectedEpsilon = 1;

    // class to sum the elements in a matrix
    private class SummationVisitor implements RealMatrixChangingVisitor
    {
        double sum = 0;

        public void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn)
        {
            sum = 0;
        }

        public double visit(int row, int column, double value) { sum += value; return value; }
        public double end() { return sum; }
    }

    /**
     * Create a UNIREP-GG test object for the specified parameters
     * @param params GLMM input parameters
     */
    public GLMMTestUnirepGeisserGreenhouse(FApproximation fMethod,
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
    public GLMMTestUnirepGeisserGreenhouse(FApproximation fMethod,
            UnivariateCdfApproximation cdfMethod, UnivariateEpsilonApproximation epsilonMethod,
            RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
            RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        super(fMethod, cdfMethod, epsilonMethod, X, XtXInverse, rank, Y,  C, U, thetaNull);
    }

    /**
     * Calculate the Geisser-Greenhouse epsilon to correct for violations
     * of sphericity
     */
    @Override
    protected void calculateEpsilon()
    {
        super.calculateEpsilon();

        if (this.epsilonMethod == UnivariateEpsilonApproximation.MULLER_BARTON_APPROX)
        {
            // calculate the expected value of the epsilon estimate
            //  E[f(lambda)] = f(lambda) + g1 / (N - r)
            //  see Muller, Barton (1989) for details
            double b = rankU;
            double g1 = 0;
            for(int i = 0; i < distinctSigmaStarEigenValues.size(); i++)
            {
                EigenValueMultiplicityPair evmI = distinctSigmaStarEigenValues.get(i);
                double firstDerivative =
                    ((2 * sumLambda)/(b * sumLambdaSquared) -
                            (2 * evmI.eigenValue * sumLambda * sumLambda)/(b*sumLambdaSquared*sumLambdaSquared));

                double secondDerivative =
                    (2 / (b * sumLambdaSquared) -
                            (8*evmI.eigenValue*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared) +
                            (8*evmI.eigenValue*evmI.eigenValue*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared*sumLambdaSquared) -
                            (2*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared));

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

            expectedEpsilon = epsilonD  + g1 / (totalN - rank);

        }
        else
        {
            // calculate the expected value of the epsilon estimate
            // see Muller, Edwards, Taylor (2004) for details

            // build a vector with each eigen value repeated per its multiplicity
            ArrayRealVector eigenColumnVector = new ArrayRealVector();
            for(EigenValueMultiplicityPair evmp: distinctSigmaStarEigenValues)
            {
                // there is probably a more memory-efficient method to do this.
                eigenColumnVector = eigenColumnVector.append(new ArrayRealVector((int) evmp.multiplicity, evmp.eigenValue));
            }

            RealMatrix outerProduct = eigenColumnVector.outerProduct(eigenColumnVector);
            double sum = outerProduct.walkInOptimizedOrder(new SummationVisitor());
            double nu = (totalN - rank);
            double expT1 = (2*nu*sumLambdaSquared) + (nu*nu*sumLambda*sumLambda);
            double expT2 = (nu*(nu + 1)*sumLambdaSquared) + (nu*sum*sum);
            expectedEpsilon = (1/(double)rankU)*(expT1/expT2);
        }
        // ensure that expected value is within bounds 1/b to 1
        if (expectedEpsilon != Double.NaN)
        {
            if (expectedEpsilon < 1/rankU)
            {
                expectedEpsilon = 1/rankU;
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
            powerNullNDFCorrection = epsilonD;
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
            powerNullDDFCorrection = expectedEpsilon;
        else
            powerNullDDFCorrection = epsilonD;

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
