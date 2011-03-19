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

import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

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
    double unirepEpsilon = Double.NaN;
    double unirepEpsilonExpectedValue = Double.NaN;

	/**
	 * Create a UNIREP-HF test object for the specified parameters
	 * @param params GLMM input parameters
	 */
    public GLMMTestUnirepHuynhFeldt(FApproximation fMethod, 
    		UnivariateCdfApproximation cdfMethod,
    		RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull, 
    		RealMatrix beta, RealMatrix sigmaError)
    {
        super(fMethod, cdfMethod, Xessence, XtXInverse, perGroupN, rank,
        		C, U, thetaNull, beta, sigmaError);

        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();

        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection();
    }
    
    /**
     * Reset the per group sample size for this test.  Recalculates epsilon
     * and expected value of epsilon
     * @param perGroupN per group sample size
     */
    @Override
    public void setPerGroupSampleSize(int perGroupN)
    {
    	super.setPerGroupSampleSize(perGroupN);
        calculateUnirepCorrection();
    }

    /**
     * Calculate the denominator degrees of freedom for the UNIREP-HF, based on
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
        
        double df = Double.NaN;

        // HF correction, we multiply the ddf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = b*(totalN - rank)*this.unirepEpsilonExpectedValue;
        else
            df = b*(totalN - rank)*this.unirepEpsilon;

        return df;
    }

    /**
     * Calculate the non-centrality parameter for the UNIREP-HF, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    @Override
    public double getNonCentrality(DistributionType type)
    {
        double a = C.getRowDimension();
        double b = U.getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*unirepEpsilon;
    }

    /**
     * Calculate the numerator degrees of freedom for the UNIREP-HF, based on
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
        
        double df = Double.NaN;

        // HF correction, we multiply the ndf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = a * b * this.unirepEpsilonExpectedValue;
        else
            df = a * b * this.unirepEpsilon;
        
        return df;
    }

    /**
     * Calculate the Huynh-Feldt epsilon to correct for violations
     * of sphericity
     */
    private void calculateUnirepCorrection()
    {          
        int b = new SingularValueDecompositionImpl(U).getRank();
        // get the sigmaStar matrix: U' *sigmaError * U
        RealMatrix sigmaStar = U.transpose().multiply(sigmaError.multiply(U));
        // ensure symmetry
        sigmaStar = sigmaStar.add(sigmaStar.transpose()).scalarMultiply(0.5); 
        // normalize
        sigmaStar = sigmaStar.scalarMultiply(1/sigmaStar.getTrace());

        // get the eigen values of the normalized sigmaStar matrix
        double[] eigenValues = new EigenDecompositionImpl(sigmaStar, TOLERANCE).getRealEigenvalues();
        if (eigenValues.length <= 0) throw new IllegalArgumentException("Failed to compute eigenvalues for sigma* matrix");
        Arrays.sort(eigenValues); // put eigenvalues in ascending order

        // calculate epsilon (correction for violation of sphericity)
        // to avoid looping over the eigenvalues twice, we also calculate the multiplicity for distinct eigenvalues

        // list of distinct eigenvalues with multiplicity
        ArrayList<EigenValueMultiplicityPair> distinctEigenValues = new ArrayList<EigenValueMultiplicityPair>();
        // initialize values for the first eigen value
        double first = eigenValues[0];
        distinctEigenValues.add(new EigenValueMultiplicityPair(first, 1));
        double sumLambda = first;
        double sumLambdaSquared = first * first;

        // loop over remaining eigen values, saving distinct eigen values
        for(int i = 1; i < eigenValues.length; i++)
        {
            double value = eigenValues[i];
            // build the sum & sum of squares of eigen values
            sumLambda += value;
            sumLambdaSquared += value * value;

            // determine if this is a distinct eigen value and calculate multiplicity
            EigenValueMultiplicityPair prev = distinctEigenValues.get(distinctEigenValues.size()-1);
            if (Math.abs(prev.eigenValue - value) > TOLERANCE)
            {
                // found new distinct eigen value
                distinctEigenValues.add(new EigenValueMultiplicityPair(value, 1));
            }
            else
            {
                // repeat of same eigenvalue, so  increment the multiplicity
                prev.multiplicity++;
            }
        }

        // calculate estimate of epsilon (correction for violation of spehericity assumption)
        unirepEpsilon = (sumLambda*sumLambda) / (b * (sumLambdaSquared));

        // For Huynh Feldt, we also need the expected value of the epsilon
        // estimate.  Note, the estimates of epsilon represent different functions of the
        // eigenvalues, so the resulting derivatives and expected values are specific to each test

        // calculate the expected value of the epsilon estimate
        // E[h(lambda)] = h(lambda) + g1 / (N - r)
        // h(lambda) = h1(lambda) / (b*h2(lambda)
        // see Muller, Barton (1989) for details
        double h1 = totalN * sumLambda * sumLambda - 2 * sumLambdaSquared;
        double h2 = (totalN - rank) * sumLambdaSquared - (sumLambda * sumLambda);
        double g1 = 0;
        for(int i = 0; i < distinctEigenValues.size(); i++)
        {
            EigenValueMultiplicityPair evmI = distinctEigenValues.get(i);
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
            for(int j = 0; j < distinctEigenValues.size(); j++)
            {
                if (i != j)
                {
                    EigenValueMultiplicityPair evmJ = distinctEigenValues.get(j);
                    // accumulate second term of g1
                    g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
                            (evmI.eigenValue - evmJ.eigenValue));
                }
            }
        }

        this.unirepEpsilonExpectedValue = (totalN*b*unirepEpsilon - 2)/(b*(totalN - rank - b*unirepEpsilon))  + g1 / (totalN - rank);


        // ensure that expected value is within bounds 1/b to 1
        if (unirepEpsilonExpectedValue != Double.NaN)
        {
            if (unirepEpsilonExpectedValue < 1/b)
            {
                unirepEpsilonExpectedValue = 1/b;
            }
            else if (unirepEpsilonExpectedValue > 1)
            {
                unirepEpsilonExpectedValue = 1;
            }
        }

    }

}
