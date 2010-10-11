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

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

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
    public GLMMTestUnirepBox(GLMMPowerParameters params)
    {
        // unirep base class will calculate epsilon for box correction
        super(params);
    }
    
    /**
     * Calculate the denominator degrees of freedom for the UNIREP-Box, based on
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

        // for the conservative, or "Box" test, we adjust by epsilon only for
        // power analysis under the alternative.  The ddf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = b*(N - r) * this.unirepEpsilon;
        else
            df = (N - r);
        
        return df;
    }

    /**
     * Calculate the non-centrality parameter for the UNIREP-Box, based on
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
     * Calculate the numerator degrees of freedom for the UNIREP-Box, based on
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

        // for the conservative, or "Box" test, we adjust by epsilon only for
        // power analysis under the alternative.  The ndf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = a * b * this.unirepEpsilon;
        else
            df = a;

        return df;
    }

}
