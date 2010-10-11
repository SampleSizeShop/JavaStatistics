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

import jsc.distributions.FishersF;

import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

/**
 * Abstract base class for statistical tests for the GLMM
 * @author Sarah Kreidler
 *
 */
public abstract class GLMMTest
{    
    // for the unirep test, the degrees of freedom change depending on two 
    // factors - data analysis (i.e. simulation or model fit) vs. power analysis,
    // and whether we need df for the F distribution under the null or the 
    // alternative hypothesis
    public enum DistributionType
    {
        POWER_NULL,
        POWER_ALTERNATIVE,
        DATA_ANALYSIS_NULL
    };
    
    // store incoming parameters
    protected GLMMPowerParameters params;

    // cache some common operations like rank, sample size, etc.
    protected double N; // total sample size (row dimension of design matrix)
    protected double r; // rank of the design matrix
    
    /**
     * Create a statistical test for the given set of GLMM parameters
     * @param params GLMM input parameters
     */
    public GLMMTest(GLMMPowerParameters params)
    {
        this.params = params;
        RealMatrix X = params.getDesign();
        N = X.getRowDimension();
        r = params.getDesignRank();
    }

    /**
     * Calculate the critical F value under the specified distribution
     * 
     * @param type distribution type
     * @param alpha type I error level
     * @return critical F
     * 
     */
    public double getCriticalF(DistributionType type, double alpha)
    throws IllegalArgumentException
    {                               
        
        double ndf = getNumeratorDF(type);
        double ddf = getDenominatorDF(type);

        FishersF centralFDist = new FishersF(ndf, ddf);
        return centralFDist.inverseCdf(1 - alpha);
    }
    
    /**
     * Prototype for getting numerator degrees of freedom for the
     * statistical test under the specified distribution type.
     * 
     * @param type distribution type
     * @return numerator degrees of freedom
     */
    abstract public double getNumeratorDF(DistributionType type);
    
    /**
     * Prototype for getting denominator degrees of freedom for the
     * statistical test under the specified distribution type.
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     */
    abstract public double getDenominatorDF(DistributionType type);
    
    /**
     * Prototype for getting observed F value for the
     * statistical test under the specified distribution type.
     * 
     * @param type distribution type
     * @return observed F
     */
    abstract public double getObservedF(DistributionType type);
    
    /**
     * Prototype for getting the non-centrality parameter for the
     * statistical test under the specified distribution type.
     * 
     * @param type distribution type
     * @return non-centrality parameter
     */
    abstract public double getNonCentrality(DistributionType type);
    
    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    protected RealMatrix getHypothesisSumOfSquares(GLMMPowerParameters params)
    {
        // convenience variables
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix B = params.getScaledBeta();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix theta0 = params.getTheta();
        
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(B.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(theta0);
        // get X'X invserse
        RealMatrix XtXInverse = params.getXtXInverse();
        
        // the middle term [C(X'X)-1C']-1
        RealMatrix cxxc = C.multiply(XtXInverse.multiply(C.transpose()));
        RealMatrix cxxcInverse = new LUDecompositionImpl(cxxc).getSolver().getInverse();
        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(cxxcInverse.multiply(thetaDiff));
        
        return hss;
        
    }
    
    /**
     * Calculate the sum of squares error matrix (the E matrix)
     * 
     * @param params matrices input by the user
     * @return error sum of squares
     */
    protected RealMatrix getErrorSumOfSquares(GLMMPowerParameters params)
    {        
        RealMatrix U = params.getWithinSubjectContrast();
        return U.transpose().multiply(params.getScaledSigmaError().multiply(U)).scalarMultiply(N-r);
    }    
    
}
