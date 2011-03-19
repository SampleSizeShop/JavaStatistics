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

/**
 * Abstract base class for statistical tests for the GLMM
 * @author Sarah Kreidler
 *
 */
public abstract class GLMMTest
{    
	// cutoff for noncentrality parameters to avoid very small negative numbers
	protected static final double TOLERANCE = 0.0000000000001;
	// type of approximation to use for unirep
	public enum UnivariateCdfApproximation
	{
		MULLER_BARTON_APPROX,
		MULLER_EDWARDS_TAYLOR_APPROX,
		MULLER_EDWARDS_TAYLOR_EXACT,
		MULLER_EDWARDS_TAYLOR_EXACT_APPROX
	};

	// available moment approximation methods
	public enum FApproximation
	{
		NONE,
		PILLAI_ONE_MOMENT,
		PILLAI_ONE_MOMENT_OMEGA_MULT,
		MCKEON_TWO_MOMENT,
		MCKEON_TWO_MOMENT_OMEGA_MULT,
		MULLER_TWO_MOMENT,
		MULLER_TWO_MOMENT_OMEGA_MULT,
		RAO_TWO_MOMENT,
		RAO_TWO_MOMENT_OMEGA_MULT
	};
	
    // for the unirep test, the degrees of freedom change depending on two 
    // factors - data analysis (i.e. simulation or model fit) vs. power analysis,
    // and whether we need df for the F distribution under the null or the 
    // alternative hypothesis
    public enum DistributionType
    {
    	DATA_ANALYSIS_NULL,
    	POWER_NULL,
    	POWER_ALTERNATIVE    	
    };
    
    // approximation information
    protected FApproximation fMethod;
    protected UnivariateCdfApproximation cdfMethod;

    // study design matrices
    protected RealMatrix Xessence;
    protected RealMatrix XtXInverse;
    protected double totalN; // total sample size 
    protected double rank; // rank of the design matrix
    // contrasts
    protected RealMatrix C; // between subject
    protected RealMatrix U; // within subject
    // null hypothesis values
    protected RealMatrix thetaNull; 
    // estimated beta parameters - only for power analysis
    protected RealMatrix beta;
    // estimated sigma error matrix
    protected RealMatrix sigmaError;
    // the M matrix: [C(X'X)-1C']-1
    protected RealMatrix M = null;
    
    /**
     * Create a statistical test for the given set of GLMM parameters
     * @param params GLMM input parameters
     */
    public GLMMTest(FApproximation fMethod, UnivariateCdfApproximation cdfMethod,
    		RealMatrix Xessence, RealMatrix XtXInverse, int perGroupN, int rank,
    		RealMatrix C, RealMatrix U, RealMatrix thetaNull, 
    		RealMatrix beta, RealMatrix sigmaError)
    {
        this.fMethod = fMethod;
        this.cdfMethod = cdfMethod;
        this.Xessence = Xessence;
        this.XtXInverse =  XtXInverse;
        this.totalN =  Xessence.getRowDimension() * perGroupN; 
        this.rank = rank; 
        this.C = C; 
        this.U =  U; 
        this.thetaNull =  thetaNull; 
        this.beta =  beta;
        this.sigmaError = sigmaError;
        
        // cache the value of M
        RealMatrix cxxcEssence = C.multiply((XtXInverse).multiply(C.transpose()));
        RealMatrix cxxcEssenceInverse = new LUDecompositionImpl(cxxcEssence).getSolver().getInverse();
        this.M = cxxcEssenceInverse.scalarMultiply(perGroupN);
    }
    
    /**
     * Reset the total sample size for this test
     * @param totalN total sample size
     */
    public void setPerGroupSampleSize(int perGroupN)
    {
    	this.totalN = Xessence.getRowDimension() * perGroupN; 
        RealMatrix cxxcEssence = C.multiply((XtXInverse).multiply(C.transpose()));
        RealMatrix cxxcEssenceInverse = new LUDecompositionImpl(cxxcEssence).getSolver().getInverse();
        this.M = cxxcEssenceInverse.scalarMultiply(perGroupN);
    }
    
    /**
     * Reset the total sample size for this test
     * @param totalN total sample size
     */
    public void setBeta(RealMatrix beta)
    {
    	this.beta = beta;
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
    protected RealMatrix getHypothesisSumOfSquares()
    {        
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(beta.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
        
        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(M.multiply(thetaDiff));
        
        return hss;
        
    }
    
    /**
     * Calculate the sum of squares error matrix (the E matrix)
     * 
     * @param params matrices input by the user
     * @return error sum of squares
     */
    protected RealMatrix getErrorSumOfSquares()
    {        
        return U.transpose().multiply(sigmaError.multiply(U)).scalarMultiply(totalN-rank);
    }    
    
}
