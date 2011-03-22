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
import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;
import org.apache.commons.math.util.MathUtils;
import edu.cudenver.bios.matrix.GramSchmidtOrthonormalization;
import edu.cudenver.bios.matrix.MatrixUtils;

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
    protected static final double TOLERANCE = 0.000000000001;
    protected double unirepEpsilon = Double.NaN;
    protected double estimatedSigmaCorrection = 1;
    protected int nuForEstimatedSigma = 0;
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
    		UnivariateCdfApproximation cdfMethod,
    		RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull, 
    		RealMatrix beta, RealMatrix sigmaError, int nuForEstimatedSigma)
    {
        super(fMethod, cdfMethod, Xessence, XtXInverse, perGroupN, rank,
        		C, U, thetaNull, beta, sigmaError);
        this.nuForEstimatedSigma = nuForEstimatedSigma;
        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();
        
        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection(nuForEstimatedSigma);
    }
    
	/**
	 * Create a UNIREP test object for data analysis
	 * @param params GLMM input parameters
	 */
    public GLMMTestUnivariateRepeatedMeasures(FApproximation fMethod, 
    		UnivariateCdfApproximation cdfMethod,
    		RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        super(fMethod, cdfMethod, X, XtXInverse, rank, Y,  C, U, thetaNull);
        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();
        
        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection(nuForEstimatedSigma);
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
        
        double df = Double.NaN;
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected test, we adjust by epsilon only for
        // power analysis under the alternative.  The ddf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = b*(totalN - rank) * this.unirepEpsilon;
        else
            df = b*(totalN - rank);
        
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
        double a = C.getRowDimension();
        double b = U.getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*(nuForEstimatedSigma > 0 ? this.estimatedSigmaCorrection : this.unirepEpsilon);
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
        
        double df = Double.NaN;
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected, we adjust by epsilon only for
        // power analysis under the alternative.  The ndf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = a * b * (nuForEstimatedSigma > 0 ? this.estimatedSigmaCorrection : this.unirepEpsilon);
        else
            df = a * b;

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
            throw new InvalidMatrixException("Failed to compute Unirep statistic: hypothesis and error matrices must be square and same dimensions");

        return H.getTrace() / E.getTrace();
    }
    
    /**
     * Calculate the epsilon to correct for violations of sphericity
     */
    private void calculateUnirepCorrection(int nuForEstimatedSigma)
    {          
    	int a = new SingularValueDecompositionImpl(C).getRank();
        int b = new SingularValueDecompositionImpl(U).getRank();
        // get the sigmaStar matrix: U' *sigmaError * U
        RealMatrix sigmaStar = U.transpose().multiply(sigmaError.multiply(U));
        double sigmaStarTrace = sigmaStar.getTrace();
        // ensure symmetry
        sigmaStar = sigmaStar.add(sigmaStar.transpose()).scalarMultiply(0.5); 
        // normalize
        sigmaStar = sigmaStar.scalarMultiply(1/sigmaStarTrace);
        
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
        
        if (nuForEstimatedSigma > 0)
        {
        	// if sigma error is estimated, we adjust slightly
        	RealMatrix estSigStar = getErrorSumOfSquares().scalarMultiply((double)1/(double)(totalN - rank));
        	double estSigStarTrace = estSigStar.getTrace();
        	double estSigStarTraceSquared = estSigStarTrace * estSigStarTrace;
        	double estSigStarSumSquares = MatrixUtils.getSumOfSquares(estSigStar);
        	// Generalized approx unbiased
        	double epsilonTilde =  ( (nuForEstimatedSigma+1) * estSigStarTraceSquared - 2 * estSigStarTraceSquared ) /
                             ( b * (nuForEstimatedSigma * estSigStarSumSquares - estSigStarTraceSquared) ) ; 
          	if (epsilonTilde > 1) epsilonTilde = 1; 
           	double mult = nuForEstimatedSigma * nuForEstimatedSigma + nuForEstimatedSigma - 2;
           	
           	RealMatrix H = getHypothesisSumOfSquares();
           	double epsilonHatNumerator =  
           		nuForEstimatedSigma * (estSigStarTraceSquared * nuForEstimatedSigma * (nuForEstimatedSigma+1) + 
           				estSigStarTrace * H.getTrace() * (2*mult/a) - 
           				estSigStarSumSquares * 2 * nuForEstimatedSigma) ;
           	double epsilonHatDenominator = (estSigStarSumSquares * nuForEstimatedSigma * nuForEstimatedSigma + 
           			(estSigStar.multiply(H)).getTrace() *(2*mult/a) - 
           			estSigStarTraceSquared* nuForEstimatedSigma)* nuForEstimatedSigma ;
           	double epsilonHat = epsilonHatNumerator / (b * epsilonHatDenominator);
           	
           	estimatedSigmaCorrection = epsilonHat;
//           	E_1_2 = EXEPS;
//            	IF (POWERCALC={6}) THEN E_1_2 = EPSTILDE_RM;  *** for HF crit val ***;
//            	IF (POWERCALC={7}) THEN E_1_2 = EPS;     *** for GG crit val ***; 
//               	IF (CDFPOWERCALC = 1) THEN E_3_5 = EPS;  *** MB ***;
//                                      ELSE E_3_5 = EPSNHAT;
//           	E_4=EPS; 
//           	CL1DF = B * NU_EST * E_4 / E_3_5;
//                END;
        	
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
        	UtU.subtract(org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(UtU.getRowDimension()));
        
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
        
        if (maxValue > MathUtils.SAFE_MIN)
        {
            // U matrix deviates from Identity, so create one that is orthonormal
            RealMatrix Utmp = new GramSchmidtOrthonormalization(U).getQ();
            U = Utmp;
        }
    }
}
