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

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;

/**
 * Implementation of the univariate approach to repeated measures test 
 * with Box correction (UNIREP-Box) for the general linear multivariate model. 
 * 
 * @see GLMMTestUnivariateRepeatedMeasures
 * @author Sarah Kreidler
 *
 */
public class GLMMTestUnirepBox extends GLMMTestUnivariateRepeatedMeasures
{
    double unirepEpsilonExpectedValue = Double.NaN;
    
	/**
	 * Create a UNIREP-Box test object for the specified parameters
	 * @param params GLMM input parameters
	 */
    public GLMMTestUnirepBox(FApproximation fMethod, 
    		UnivariateCdfApproximation cdfMethod,
    		RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull, 
    		RealMatrix beta, RealMatrix sigmaError, int nuEst)
    {
        // unirep base class will calculate epsilon for box correction
        super(fMethod, cdfMethod, Xessence, XtXInverse, perGroupN, rank,
        		C, U, thetaNull, beta, sigmaError, nuEst);
    }
    
	/**
	 * Create a UNIREP test object for data analysis
	 * @param params GLMM input parameters
	 */
    public GLMMTestUnirepBox(FApproximation fMethod, 
    		UnivariateCdfApproximation cdfMethod,
    		RealMatrix X, RealMatrix XtXInverse, int rank, RealMatrix Y,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull)
    {
        super(fMethod, cdfMethod, X, XtXInverse, rank, Y,  C, U, thetaNull);
    }
    
    /**
     * Calculate the correction factors for numerator degrees of 
     * freedom for data analysis, power under the null and power
     * under the alternative
     */
    @Override
    protected void calculateNDFCorrection()
    {
    	//double b = (double) U.getColumnDimension();
        dataAnalysisNDFCorrection = 1.0/rankU;
        powerNullNDFCorrection = 1.0/rankU;
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
    @Override
    protected void calculateDDFCorrection()
    {
        dataAnalysisDDFCorrection = 1.0/rankU;
        powerNullDDFCorrection = 1.0/rankU;
        powerAlternativeDDFCorrection = epsilonD;
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
