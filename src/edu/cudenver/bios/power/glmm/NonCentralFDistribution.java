/*
* GLIMMPSE (General Linear Multivariate Model Power and Sample size Estimator)
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
import jsc.distributions.NoncentralFishersF;
import jsc.distributions.Normal;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;

/**
 * 
 * This class calculates Pr(Fcrit < F(ndf,ddf,nonCentrality)) using one of four methods. 
 * The first, most common method uses the cumulative distribution function of the non-central F.  
 * If the CDF method will fail, the second method uses the Tiku approximation to the non-central F 
 * (see Johnson and Kotz).  In situations where the TIKU will fail or be inaccurate, a 
 * normal approximation is used.
 *  
 * @author translated from POWERLIB (SAS/IML) written by Keith Muller
 *
 */
public class NonCentralFDistribution
{     
    private static final int STARTING_F = 10;
    
    // supported approximation methods for the F distribution
    protected enum FMethod
    {
        CDF,  // use the regular non-central F distribution
        TIKU_APPROXIMATION, // use the Tiku approximation to the non-central F
        NORMAL_APPROXIMATION // approximate with a normal distribution
    };
    
    // numerator and denominator degrees of freedom
    protected double ndf;
    protected double ddf;
    protected double nonCentrality;
    // method of approximation
    protected FMethod method = FMethod.CDF;
    
    // non-central F distribution for use with the CDF method
    protected NoncentralFishersF nonCentralF= null;
    
    // Intermediate parameters for use with the Tiku approximation method
    protected double TikuC;
    protected double TikuB;
    
    // normal distribution for use with the NORMAL_APPROXIMATION method
    protected Normal normal = null;

    /**
     * Function calculating the difference between the probability of a target quantile 
     * and the  (used by the bisection solver from Apache Commons Math)
     * @see org.apache.commons.math.analysis.UnivariateRealFunction
     */
    private class NonCentralFQuantileFunction implements UnivariateRealFunction
    {
        protected double quantile;
        
        public NonCentralFQuantileFunction(double quantile)
        {
            this.quantile = quantile;
        }
        
        public double value(double n)
        {
            return cdf(n) - quantile;
        }
    }
    
    /**
     * Create a non-central F distribution with the specified numerator degrees of freedom,
     * denominator degrees of freedom, and non-centrality parameter
     * 
     * @param numeratorDegreesOfFreedom
     * @param denominatorDegreesOfFreedom
     * @param nonCentralityParameter
     */
    public NonCentralFDistribution(double numeratorDegreesOfFreedom, 
            double denominatorDegreesOfFreedom, double nonCentralityParameter)
    {
        ndf = numeratorDegreesOfFreedom;
        ddf = denominatorDegreesOfFreedom;
        nonCentrality = nonCentralityParameter;
        
        // set the method of approximation and initialize appropriately
        if ((ndf < Math.pow(10, 4.4) && ddf < Math.pow(10, 5.4) && 
                this.nonCentrality < Math.pow(10, 6.4)) ||
                (ndf < Math.pow(10, 6) && ddf < 10 && nonCentrality <= Math.pow(10, 6)))
        {
            // cdf valid in these ranges
            method = FMethod.CDF;
            nonCentralF = new NoncentralFishersF(ndf, ddf, nonCentrality);
        }
        else if (ndf < Math.pow(10, 9.2) && ndf >= 1 &&
                ddf < Math.pow(10, 9.2) && ddf > Math.pow(10, 0.6) &&
                this.nonCentrality < Math.pow(10, 6.4))
        {
            // use Tiku approximation for extreme ndf values
            method = FMethod.TIKU_APPROXIMATION;
            double TikuH = 2 * Math.pow(ndf + nonCentrality, 3) +  
                3 * (ndf + nonCentrality) * (ndf + 2 * nonCentrality) *(ddf - 2) + 
                (ndf + 3 * nonCentrality) * Math.pow(ddf - 2, 2);
            double TikuK = Math.pow(ndf + nonCentrality, 2) + 
                (ddf - 2) * (ndf + 2 * nonCentrality);
            double TikuNdf = Math.floor(0.5 * (ddf - 2)*(Math.sqrt((TikuH*TikuH) / 
                    ((TikuH*TikuH) - 4 * Math.pow(TikuK, 3))) - 1));
            TikuC = (TikuNdf / ndf) / (2 * TikuNdf + ddf - 2) * (TikuH / TikuK);
            TikuB = - ddf / (ddf - 2) * (TikuC - 1 - nonCentrality / ndf);
            nonCentralF = new NoncentralFishersF(TikuNdf, ddf, 0);
        }
        else
        {
            method = FMethod.NORMAL_APPROXIMATION;
            normal = new Normal();
        }
    }

    /**
     * For this non-central F distribution, F, this method returns P(F < f). 
     * 
     * @param Fcritical critical point for which to calculate cumulative probability
     * @return
     */
    public double cdf(double Fcritical)
    {
        if (Fcritical <= 0) return 0;
        switch (method)
        {
        case TIKU_APPROXIMATION:
            double TikuFcrit = (Fcritical - TikuB) / TikuC;
            if (TikuFcrit <= 0) return 0;
            return nonCentralF.cdf(TikuFcrit);
        case NORMAL_APPROXIMATION:
            double p1 = 1 / 3;
            double p2 = -2;
            double p3 = 1 / 2;
            double p4 = 2 / 3;
            double arg1 = ((ndf * Fcritical) / (ndf + nonCentrality));
            double arg2 = (2 / 9) * (ndf + (2 * nonCentrality)) * (Math.pow(ndf + nonCentrality, p2));
            double arg3 = (2 / 9) * (1 / ddf);
            double numZ = Math.pow(arg1, p1) - (arg3 * Math.pow(arg1, p1)) - (1 - arg2);
            double denZ =  Math.pow((arg2 + arg3 * Math.pow(arg1, p4)), p3);
            double zScore = numZ / denZ;
            double absZScore = Math.abs(zScore);
            if (absZScore < 6)
            {
                return normal.cdf(zScore);
            }
            else
            {
                if (zScore < -6) 
                    return 0;
                else
                    return 1; // Z > 6
            }
        default: // method = CDF
            return nonCentralF.cdf(Fcritical);
        }      

    }
    
    /**
     * For this non-central F distribution, F, this function returns the critical value, f,
     * such that P(F < f). 
     * 
     * @param quantile quantile
     * @return
     */
    public double inverseCDF(double quantile)
    {
        if (quantile <= 0) return Double.NaN;
        if (quantile >= 1) return Double.POSITIVE_INFINITY;
        
        if (method == FMethod.CDF && nonCentrality == 0)
        {
            // inverseCdf throws an illegal argument exception when the
            // non-centrality parameter is non-zero.  So we just bisection solve 
            // unless we're really dealing with a central F.  Ah, the hazards of
            // pulling code off of the interwebs.
            return nonCentralF.inverseCdf(quantile);
        }
        else
        {
            UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
            UnivariateRealSolver solver = factory.newBisectionSolver();

            NonCentralFQuantileFunction quantFunc = new NonCentralFQuantileFunction(quantile);

            int upperBound = findUpperBound(quantile);
            try
            {
                return solver.solve(quantFunc, 0, upperBound);
            }
            catch (Exception e)
            {
                throw new IllegalArgumentException("Failed to determine F quantile: " + e.getMessage());
            }
        }
    }
        
    private int findUpperBound(double quantile)
    {
        int upperBound = STARTING_F;

        for(double currentQuantile = 0.0; currentQuantile < quantile; upperBound *= 2)
        {
            currentQuantile = cdf(upperBound);
        }
        return upperBound;
    }
}
