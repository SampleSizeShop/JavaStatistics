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

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealMatrixChangingVisitor;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

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
    double unirepEpsilon = Double.NaN;
    double unirepEpsilonExpectedValue = Double.NaN;

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
    public GLMMTestUnirepGeisserGreenhouse(GLMMPowerParameters params)
    {
        super(params);

        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();

        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection();
    }

    /**
     * Calculate the denominator degrees of freedom for the UNIREP-GG, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix U = params.getWithinSubjectContrast();
        
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();

        double df = Double.NaN;
        
        // GG correction, we multiply the ddf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = b*(N - r)*this.unirepEpsilonExpectedValue;
        else
            df = b*(N - r)*this.unirepEpsilon; 
        
        return df;
    }

    /**
     * Calculate the non-centrality parameter for the UNIREP-GG, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    @Override
    public double getNonCentrality(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*unirepEpsilon;
    }

    /**
     * Calculate the numerator degrees of freedom for the UNIREP-GG, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        double df = Double.NaN;

        // GG correction, we multiply the ndf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = a * b * this.unirepEpsilonExpectedValue;
        else
            df = a * b * this.unirepEpsilon;
        
        return df;
    }

    /**
     * Calculate the Geisser-Greenhouse epsilon to correct for violations
     * of sphericity
     */
    private void calculateUnirepCorrection()
    {          
        RealMatrix U = params.getWithinSubjectContrast();
        int b = new SingularValueDecompositionImpl(U).getRank();
        // get the sigmaStar matrix: U' *sigmaError * U
        RealMatrix sigmaStar = U.transpose().multiply(params.getScaledSigmaError().multiply(U));
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

//         calculate the expected value of the epsilon estimate
//         E[f(lambda)] = f(lambda) + g1 / (N - r)
//         see Muller, Barton (1989) for details
//        double g1 = 0;
//        for(int i = 0; i < distinctEigenValues.size(); i++)
//        {
//            EigenValueMultiplicityPair evmI = distinctEigenValues.get(i);
//            double firstDerivative = 
//                ((2 * sumLambda)/(b * sumLambdaSquared) - 
//                        (2 * evmI.eigenValue * sumLambda * sumLambda)/(b*sumLambdaSquared*sumLambdaSquared));
//
//            double secondDerivative = 
//                (2 / (b * sumLambdaSquared) - 
//                        (8*evmI.eigenValue*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared) + 
//                        (8*evmI.eigenValue*evmI.eigenValue*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared*sumLambdaSquared) - 
//                        (2*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared)); 
//
//            // accumulate the first term of g1 (sum over distinct eigen vals of 1st derivative * eigen val ^2 * multiplicity)
//            g1 += secondDerivative * evmI.eigenValue * evmI.eigenValue * evmI.multiplicity;
//            // loop over elements not equal to current
//            for(int j = 0; j < distinctEigenValues.size(); j++)
//            {
//                if (i != j)
//                {
//                    EigenValueMultiplicityPair evmJ = distinctEigenValues.get(j);
//                    // accumulate second term of g1
//                    g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
//                            (evmI.eigenValue - evmJ.eigenValue));
//                }
//            }
//        }
//
//        this.unirepEpsilonExpectedValue = unirepEpsilon  + g1 / (N - r);

        // calculate the expected value of the epsilon estimate
        // see Muller, Edwards, Taylor (2004) for details
        
        // build a vector with each eigen value repeated per its multiplicity
        ArrayRealVector eigenColumnVector = new ArrayRealVector();
        for(EigenValueMultiplicityPair evmp: distinctEigenValues) 
        {
            // there is probably a more memory-efficient method to do this. 
            eigenColumnVector = eigenColumnVector.append(new ArrayRealVector((int) evmp.multiplicity, evmp.eigenValue));
        }

        RealMatrix outerProduct = eigenColumnVector.outerProduct(eigenColumnVector);
        double sum = outerProduct.walkInOptimizedOrder(new SummationVisitor());
        double nu = (N - r);
        double expT1 = (2*nu*sumLambdaSquared) + (nu*nu*sumLambda*sumLambda);
        double expT2 = (nu*(nu + 1)*sumLambdaSquared) + (nu*sum*sum);
        unirepEpsilonExpectedValue = (1/(double)b)*(expT1/expT2);
  
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
