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

import jsc.distributions.FishersF;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.distribution.ChiSquareTerm;
import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.distribution.WeightedSumOfNoncentralChiSquaresDistribution;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

/**
 * Class representing the distribution of the non-centrality parameter in
 * the general linear multivariate model.  Used by the GLMMPowerCalculator class 
 * for computing unconditional and quantile power.
 * 
 * @see edu.cudenver.bios.power.GLMMPowerCalculator
 * @author Sarah Kreidler
 */
public class NonCentralityDistribution
{
    private static final double TOLERANCE = 0.000000000001;
    private static final double ACCURACY = 0.001;
    // intermediate forms 
    protected RealMatrix T1 = null;
    protected RealMatrix FT1 = null;
    protected RealMatrix S = null;
    protected RealMatrix mzSq = null;
    protected double H1;
    protected double H0 = 0;
    int qF;
    int a;
    int N;
    double[] sEigenValues;
    int sStar = 0;
    // indicates if an "exact" cdf should be calculated via Davie's algorithm or
    // with the Satterthwaite approximation from Glueck & Muller
    protected boolean exact;

    /**
     * Function calculating the difference between the probability of a target quantile 
     * and the  (used by the bisection solver from Apache Commons Math)
     * @see org.apache.commons.math.analysis.UnivariateRealFunction
     */
    private class NonCentralityQuantileFunction implements UnivariateRealFunction
    {
        protected double quantile;
        
        public NonCentralityQuantileFunction(double quantile)
        {
            this.quantile = quantile;
        }
        
        public double value(double n)
        {
            return cdf(n) - quantile;
        }
    }
    
    /**
     * Create a non-centrality distribution for the specified inputs.
     * @param params GLMM input parameters
     * @param exact if true, Davie's algorithm will be used to compute the cdf, 
     * otherwise a Satterthwaite style approximation is used.
     * @throws IllegalArgumentException
     */
    public NonCentralityDistribution(Test test, RealMatrix F, RealMatrix FtFinverse, int N, 
    		FixedRandomMatrix CFixedRand, RealMatrix U, 
    		RealMatrix thetaNull, RealMatrix beta, 
    		RealMatrix sigmaError, RealMatrix sigmaG, boolean exact)
    throws IllegalArgumentException
    {
        this.exact = exact;
        try
        {                    
            // TODO: need to calculate H0, need to adjust H1 for Unirep
            // get design matrix for fixed parameters only
//            RealMatrix F = params.getDesignEssence().getFullDesignMatrixFixed();
            qF = F.getColumnDimension();
            a = CFixedRand.getCombinedMatrix().getRowDimension();
//            N = F.getRowDimension();
            // get fixed contrasts
            RealMatrix Cfixed = CFixedRand.getFixedMatrix();
            RealMatrix CGaussian = CFixedRand.getRandomMatrix();
            // build intermediate terms h1, S
            if (FtFinverse == null)
            {
            	FtFinverse = new LUDecompositionImpl(F.transpose().multiply(F)).getSolver().getInverse();
            }
            RealMatrix P = Cfixed.multiply(FtFinverse).multiply(F.transpose());
            RealMatrix PPt = P.multiply(P.transpose());            
            T1 = new LUDecompositionImpl(PPt).getSolver().getInverse();
            FT1 = new CholeskyDecompositionImpl(T1).getL();
            // calculate theta difference
//            RealMatrix theta0 = params.getTheta();
            RealMatrix C = CFixedRand.getCombinedMatrix();
//            RealMatrix B = params.getScaledBeta();
//            RealMatrix U = params.getWithinSubjectContrast();
            // thetaHat = C * Beta * U
            RealMatrix thetaHat = C.multiply(beta.multiply(U));
            // thetaHat - thetaNull.  
            RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
            // TODO: specific to HLT or UNIREP
            RealMatrix sigmaStarInverse = getSigmaStarInverse(U, sigmaError, test);
            RealMatrix H1matrix = thetaDiff.transpose().multiply(T1).multiply(thetaDiff).multiply(sigmaStarInverse);
            H1 = H1matrix.getTrace();
            // matrix which represents the non-centrality parameter as a linear combination of chi-squared r.v.'s
            S = FT1.transpose().multiply(thetaDiff).multiply(sigmaStarInverse).multiply(thetaDiff.transpose()).multiply(FT1).scalarMultiply(1/H1);
            // we use the S matrix to generate the F-critical, numerical df's, and denominator df's 
            // for a central F distribution.  The resulting F distribution is used as an approximation
            // for the distribution of the non-centrality parameter
            // See formulas 18-21 and A8,A10 from Glueck & Muller (2003) for details
            EigenDecompositionImpl sEigenDecomp = new EigenDecompositionImpl(S, TOLERANCE);
            sEigenValues = sEigenDecomp.getRealEigenvalues();
            // calculate H0
            if (sEigenValues.length > 0) H0 = H1 * (1 - sEigenValues[0]);
            if (H0 <= 0) H0 = 0;
            
            // count the # of positive eigen values
            for(double value: sEigenValues) 
            {
                if (value > 0) sStar++;
            }
            // TODO: throw error if sStar is <= 0
            double stddevG = Math.sqrt(sigmaG.getEntry(0, 0));
            RealMatrix svec = sEigenDecomp.getVT();
            mzSq = svec.multiply(FT1.transpose()).multiply(CGaussian).scalarMultiply(1/stddevG);
            for(int i = 0; i < mzSq.getRowDimension(); i++)
            {
                for (int j = 0; j < mzSq.getColumnDimension(); j++)
                {
                    double entry = mzSq.getEntry(i, j);
                    mzSq.setEntry(i, j, entry*entry); // TODO: is there an apache function to do this?
                }
            }
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e.getMessage());
        }
    }
    
    /**
     * Calculate the probability P(W < w), where W follows the distribution of
     * the non-centrality parameter 
     * 
     * @param w critical point for which to calculate cumulative probability
     * @return P(W < w)
     */
    public double cdf(double w)
    {
        if (H1 <= 0 || w <= H0) return 0;
        if (H1 - w <= 0) return 1;
        ArrayList<ChiSquareTerm> chiSquareTerms = new ArrayList<ChiSquareTerm>();
        
        try
        {
            double b0 = 1 - w / H1;
            double m1Positive = 0;
            double m1Negative = 0;
            double m2Positive = 0;
            double m2Negative = 0;
            //
            int numPositive = 0;
            int numNegative = 0;
            double nu;
            double delta;
            double lambda;
            double lastPositiveNoncentrality = 0; // for special cases
            double lastNegativeNoncentrality = 0; // for special cases
            
            // add in the first chi-squared term in the estimate of the non-centrality
            // (expressed as a sum of weighted chi-squared r.v.s)
            // initial chi-square term is central (delta=0) with N-qf df, and lambda = b0
            nu = N - qF;
            lambda = b0;
            delta = 0;
            chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
            // accumulate terms
            if (lambda > 0)
            {
                // positive terms
                numPositive++;
                lastPositiveNoncentrality = delta;
                m1Positive += lambda * (nu + delta);
                m2Positive += lambda * lambda * 2* (nu + 2*delta);
            }
            else if (lambda < 0)
            {
                // negative terms - we take absolute value of lambda where needed
                numNegative++;
                lastNegativeNoncentrality = delta;
                m1Negative += -1 * lambda * (nu + delta);
                m2Negative += lambda * lambda * 2* (nu + 2*delta);
            }
            
            // accumulate the remaining terms 
            for(int k = 0; k < sStar; k++)
            {
                if (k < sStar)
                {
                    // for k = 1 (well, 0 in java array terms and 1 in the paper) to sStar, chi-square term is 
                    // non-central (delta = mz^2), 1 df, lambda = (b0 - kth eigen value of S)
                    nu = 1;
                    lambda = b0 - sEigenValues[k];
                    delta = mzSq.getEntry(k, 0);
                    chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
                }
                else
                {
                    // for k = sStar+1 to a, chi-sqaure term is non-central (delta = mz^2), 1 df,
                    // lambda = b0
                    nu = 1;
                    lambda = b0;
                    delta = mzSq.getEntry(k, 0);
                    chiSquareTerms.add(new ChiSquareTerm(lambda, nu, delta));
                }
                // accumulate terms
                if (lambda > 0)
                {
                    // positive terms
                    numPositive++;
                    lastPositiveNoncentrality = delta;
                    m1Positive += lambda * (nu + delta);
                    m2Positive += lambda * lambda * 2* (nu + 2*delta);
                }
                else if (lambda < 0)
                {
                    // negative terms - we take absolute value of lambda where needed
                    numNegative++;
                    lastNegativeNoncentrality = delta;
                    m1Negative += -1 * lambda * (nu + delta);
                    m2Negative += lambda * lambda * 2* (nu + 2*delta);
                }
                // Note, we deliberately ignore terms for which lambda == 0
            }

            // handle special cases
            if (numNegative == 0) return 0;
            if (numPositive == 0) return 1;
            
            // special cases
            if (numNegative == 1 && numPositive == 1)
            {
            	double Nstar = N - qF + a - 1;
            	double Fstar = w / (Nstar * (H1 - w));
            	if (lastPositiveNoncentrality >= 0 && lastNegativeNoncentrality == 0)
            	{
                	// handle special case: CGaussian = 0, s* = 1
                    NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(Nstar, 1, lastPositiveNoncentrality);
            		return nonCentralFDist.cdf(Fstar);
            	}
            	else if (lastPositiveNoncentrality == 0 && lastNegativeNoncentrality > 0)
            	{
                	// handle special case: CGaussian = 1
                    NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(1,Nstar, lastNegativeNoncentrality);
            		return 1 - nonCentralFDist.cdf(1/Fstar);
            	}
            }
            
            if (exact)
            {
            	WeightedSumOfNoncentralChiSquaresDistribution dist	= 
            		new WeightedSumOfNoncentralChiSquaresDistribution(chiSquareTerms, ACCURACY);
            	return dist.cdf(0);
            }
            else
            {
            	// handle general case - Satterthwaite approximation
            	double nuStarPositive = 2 * (m1Positive * m1Positive) / m2Positive;
            	double nuStarNegative = 2 * (m1Negative * m1Negative) / m2Negative;
            	double lambdaStarPositive = m2Positive / (2 * m1Positive);
            	double lambdaStarNegative =  m2Negative / (2 * m1Negative);

            	// create a central F to approximate the distribution of the non-centrality parameter
            	FishersF centralFDist = new FishersF(nuStarPositive, nuStarNegative);
            	// return power based on the non-central F
            	return centralFDist.cdf((nuStarNegative*lambdaStarNegative)/(nuStarPositive*lambdaStarPositive));
            }
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
    }
    
    /**
     * For this non-centrality distribution, W, this function returns the critical value, w,
     * such that P(W < w). 
     * 
     * @param probability desired value of P(W < w)
     * @return critical w such that P(W < w)
     */
    public double inverseCDF(double probability)
    {
        if (H1 <= 0) return 0;

        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        NonCentralityQuantileFunction quantFunc = new NonCentralityQuantileFunction(probability);
        
        try
        {
            return solver.solve(quantFunc, H0, H1);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to determine non-centrality quantile: " + e.getMessage());
        }
    }
    
    /**
     * Calculate the inverse of the sigma star matrix
     * @param params GLMM input parameters
     * @return sigma star inverse
     */
    private RealMatrix getSigmaStarInverse(RealMatrix U, RealMatrix sigmaError, Test test)
    {
        // sigma* = U'*sigmaE*U
        RealMatrix sigmaStar = U.transpose().multiply(sigmaError).multiply(U);
        
        if (test == Test.HOTELLING_LAWLEY_TRACE)
        {
            return new LUDecompositionImpl(sigmaStar).getSolver().getInverse();
        }
        else
        {
            // stat should only be UNIREP (uncorrected, box, GG, or HF) at this point  
        	// (exception is thrown by valdiateParams otherwise)
            int b = sigmaStar.getColumnDimension();
            // get discrepancy from sphericity for unirep test
            double sigmaStarTrace = sigmaStar.getTrace();
            double sigmaStarSquaredTrace = sigmaStar.multiply(sigmaStar).getTrace();
            double epsilon = (sigmaStarTrace*sigmaStarTrace) / ((double) b * sigmaStarSquaredTrace);
            RealMatrix identity = MatrixUtils.createRealIdentityMatrix(b);
            return identity.scalarMultiply((double) b * epsilon / sigmaStarTrace);
        }
    }

    /**
     * Get the upper intergral bound
     * @return H1
     */
    public double getH1()
    {
        return H1;
    }   
    
    /**
     * Get the lower intergral bound
     * @return H0
     */
    public double getH0()
    {
        return H0;
    }   
}
