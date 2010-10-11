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

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.MomentApproximationMethod;

/**
 * Implementation of the Hotelling Lawley Trace (HLT) test for the
 * general linear multivariate model
 * 
 * @author Sarah Kreidler
 *
 */
public class GLMMTestHotellingLawley extends GLMMTest
{    
	/**
	 * Create a Hotelling-Lawley Trace test object for the specified parameters
	 * @param params GLMM input parameters
	 */
    public GLMMTestHotellingLawley(GLMMPowerParameters params)
    {
        super(params);
    }
    
    /**
     * Calculate the denominator degrees of freedom for the HLT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        double a = (double) C.getRowDimension();
        // b = #columns in within subject contrast matrix
        double b = (double) U.getColumnDimension();
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * ((N - r) - b -1) + 2;
        }
        else
        {
            double t1 = (N - r) * (N - r) - (N - r) * (2 * b + 3) + b * (b + 3);
            double t2 = ((N - r) * (a  + b + 1) - (a + 2 * b + b * b - 1));
            df = 4 + (a * b + 2) * (t1/t2);
        }
        // TODO Auto-generated method stub
        return df;
    }

    /**
     * Calculate the non-centrality parameter for the HLT, based on
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
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(params);
        RealMatrix errorSumOfSquares = getErrorSumOfSquares(params);
        
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix B = params.getScaledBeta();
        
        // check if we are uni or multi variate
        double p = B.getColumnDimension();
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
        
        if ((s == 1 && p > 1) ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                params.getMomentMethod() == MomentApproximationMethod.MCKEON_TWO_MOMENT_OMEGA_MULT)
        {
            HLT *= ((double)(N - r)/(double)N);
            return N * s * HLT / s;
        }
        else
        {
            return getDenominatorDF(type) * HLT / s;
        }
    }

    /**
     * Calculate the numerator degrees of freedom for the HLT, based on
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
        
        return a * b;
    }

    /**
     * Calculate the observed F for the HLT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return observed F
     * @throws IllegalArgumentException
     */
    @Override
    public double getObservedF(DistributionType type)
    {       
        if (type == DistributionType.DATA_ANALYSIS_NULL)
        {
            RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(params);
            RealMatrix errorSumOfSquares = getErrorSumOfSquares(params);
            double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
            double ddf = getDenominatorDF(type);
            double ndf = getNumeratorDF(type);
            double b = params.getWithinSubjectContrast().getColumnDimension();
            return HLT * (((N-r)-b-1)*ddf) / (ndf*(ddf-2));
        }
        else
        {
            return getNonCentrality(type) / getNumeratorDF(type);
        }
    }

    /**
     * Compute a Hotelling-Lawley Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getHotellingLawleyTrace(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Hotelling Lawley Trace: hypothesis and error matrices must be square and same dimensions");
        
        RealMatrix inverseE = new LUDecompositionImpl(E).getSolver().getInverse();
        RealMatrix HinverseE = H.multiply(inverseE);
        
        return HinverseE.getTrace();
    }
}
