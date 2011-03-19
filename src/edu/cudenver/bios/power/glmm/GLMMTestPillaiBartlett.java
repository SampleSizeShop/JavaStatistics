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

import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

/**
 * Implementation of the Pillai Bartlett Trace (PBT) test for the
 * general linear multivariate model
 * 
 * @author Sarah Kreidler
 *
 */
public class GLMMTestPillaiBartlett extends GLMMTest
{
	/**
	 * Create a Pillai Bartlett Trace test object for the specified parameters
	 * @param params GLMM input parameters
	 */
    public GLMMTestPillaiBartlett(FApproximation fMethod,
    		RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull, 
    		RealMatrix beta, RealMatrix sigmaError)
    {
        super(fMethod, null, Xessence, XtXInverse, perGroupN, rank,
        		C, U, thetaNull, beta, sigmaError);
    }   
    
    /**
     * Calculate the denominator degrees of freedom for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getDenominatorDF(DistributionType type)
    {        
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix
        double b = U.getColumnDimension();
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (fMethod == FApproximation.PILLAI_ONE_MOMENT ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * ((totalN - rank) - b + s);
        }
        else
        {
            double mu1= a * b / (totalN - rank + a);
            double factor1 = (totalN - rank + a - b) / (totalN - rank + a - 1);
            double factor2 = (totalN - rank) / (totalN - rank + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((totalN - rank + a)*(totalN - rank + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * (m1 - m2) * (1 - m1) / denom;
        }

        return df;
    }

    /**
     * Calculate the non-centrality parameter for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    @Override
    public double getNonCentrality(DistributionType type)
    {
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares();
        RealMatrix errorSumOfSquares = getErrorSumOfSquares();
        
        // check if we are uni or multi variate
        double p = beta.getColumnDimension();
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double PB = getPillaiBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
        
        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                fMethod == FApproximation.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            return totalN * s * PB / (s - PB);
        }
        else
        {
            return getDenominatorDF(type) * PB / (s - PB);
        }
    }

    /**
     * Calculate the numerator degrees of freedom for the PBT, based on
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
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (
        		fMethod == FApproximation.PILLAI_ONE_MOMENT ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = a * b;
        }
        else
        {
            double mu1= a * b / (totalN - rank + a);
            double factor1 = (totalN - rank + a - b) / (totalN - rank + a - 1);
            double factor2 = (totalN - rank) / (totalN - rank + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((totalN - rank + a)*(totalN - rank + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * m1 * (m1 - m2) / denom;
        }
        
        return df;
    }

    /**
     * Calculate the observed F for the PBT, based on
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
        
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double association = 0.0;
        
        double PB = getPillaiBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
        association = PB / s;
        double ddf = getDenominatorDF(type);
        double ndf = getNumeratorDF(type);

        return ((association) / ndf) / ((1 - association) / ddf);
    }
    
    /**
     * Compute a Pillai Bartlett Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getPillaiBartlettTrace(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Pillai-Bartlett Trace: hypothesis and error matrices must be square and same dimensions");
        
        double a = C.getRowDimension();
        double b = U.getColumnDimension();
        double s = (a < b) ? a : b;  
        double p = beta.getColumnDimension();
        
        RealMatrix adjustedH = H;
        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                fMethod == FApproximation.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            adjustedH = H.scalarMultiply(((double)(totalN - rank)/(double)totalN));
        }
            
        RealMatrix T = adjustedH.add(E);
        RealMatrix inverseT = new LUDecompositionImpl(T).getSolver().getInverse();

        RealMatrix HinverseT = adjustedH.multiply(inverseT);
        
        return HinverseT.getTrace();
    }
}
