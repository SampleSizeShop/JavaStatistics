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
package edu.cudenver.bios.power;

import java.util.ArrayList;
import java.util.List;

import jsc.distributions.FishersF;
import jsc.distributions.Normal;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.Test;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.NonCentralFDistribution;
import edu.cudenver.bios.power.glmm.NonCentralityDistribution;

/**
 * Power calculator implementation for the general linear multivariate model
 * 
 * @see PowerCalculator
 * @author Sarah Kreidler
 *
 */
public class GLMMPowerCalculator implements PowerCalculator
{
	private static final double RELATIVE_SYMMETRY_THRESHOLD = 
	    CholeskyDecompositionImpl.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;
	private static final double POSITIVITY_THRESHOLD = 
	    CholeskyDecompositionImpl.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
	    
    private static final int STARTING_SAMPLE_SIZE = 1000;
    private static final int STARTING_DETECTABLE_DIFFERENCE = 100;
    private static final int SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL = 1000;
    /**
     * Simple class to contain sample size and actual power
     * from a sample size calculation
     */
    private class SampleSize
    {
        public int sampleSize; // total sample size
        public double actualPower; // calculated power associated with total sample size
        
        public SampleSize(int sampleSize, double actualPower)
        {
            this.sampleSize = sampleSize;
            this.actualPower = actualPower;
        }
    }
    
    /**
     * Function used with Apache's bisection solver to determine the 
     * per-group sample size which most closely achieves the desired power
     */
    private class SampleSizeFunction implements UnivariateRealFunction
    {
    	GLMMPowerParameters params;

        public SampleSizeFunction(GLMMPowerParameters params)
        {
            this.params = params;
        }
        
        public double value(double n)
        {
            try
            {
                params.setGroupSampleSize((int) n); 
                double calculatedPower = getPowerByType(params);
                return params.getCurrentPower() - calculatedPower;
            }
            catch (Exception e)
            {   
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a negative number
                return Double.NaN;  
            }
        }
    }
    
    /**
     * Simple class to contain sample size and actual power
     * from a sample size calculation
     */
    private class DetectableDifference
    {
        public double betaScale; // scale factor for beta matrix
        public double actualPower; // calculated power associated with total sample size
        
        public DetectableDifference(double betaScale, double actualPower)
        {
            this.betaScale = betaScale;
            this.actualPower = actualPower;
        }
    }
    
    /**
     * Function used with Apache's bisection solver to determine the 
     * per-group sample size which most closely achieves the desired power
     */
    private class DetectableDifferenceFunction implements UnivariateRealFunction
    {
        GLMMPowerParameters params;
        
        public DetectableDifferenceFunction(GLMMPowerParameters params)
        {
            this.params = params;
        }
        
        public double value(double betaScale)
        {
            try
            {
                params.scaleBeta(betaScale);
                double calculatedPower = getPowerByType(params);
                return params.getCurrentPower() - calculatedPower;
            }
            catch (Exception e)
            {   
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a negative number
                return Double.NaN;  
            }
        }
    }
    /**
     * Class passed into Apache's TrapezoidIntegrator function to compute
     * unconditional power
     */
    private class UnconditionalPowerIntegrand implements UnivariateRealFunction
    {
        protected NonCentralityDistribution nonCentralityDist;
        protected double Fcrit;
        protected double ndf;
        protected double ddf;

        public UnconditionalPowerIntegrand(NonCentralityDistribution nonCentralityDist,
                double Fcrit, double ndf, double ddf)
        {
            this.nonCentralityDist = nonCentralityDist;
            this.Fcrit = Fcrit;
            this.ndf = ndf;
            this.ddf = ddf;
        }
        
        public double value(double t)
        {
            NonCentralFDistribution FdistTerm1 = new NonCentralFDistribution(ndf, ddf, t);
            NonCentralFDistribution FdistTerm2 = new NonCentralFDistribution(ndf+2, ddf, t);

            return nonCentralityDist.cdf(t)*(FdistTerm1.cdf(Fcrit) - FdistTerm2.cdf((Fcrit*ndf)/(ndf+2)));
        }
    }
	
	/********* public methods for the power API ************/
	
	/**
	 * Calculate a list of power values using the methodology described in
	 * Glueck & Muller 2003.
	 * If the parameters contain lists of possible scale factors, statistical
	 * tests, etc., then a power will be returned for each combination
	 * of these factors.
	 * 
	 * @see GLMMPowerParameters
	 * @see GLMMPower
	 * @param powerParams inputs to the power calculation
	 * @return list of calculated power values.
	 */
	@Override
	public List<Power> getPower(PowerParameters powerParams)
	{
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);
        
        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();
                
        // calculate the power for all variations of the study design
        for(GLMMPowerParameters.Test test = params.getFirstTest(); test != null;
        test = params.getNextTest())
        {            
        	for(GLMMPowerParameters.PowerMethod method = params.getFirstPowerMethod(); method != null;
        	method = params.getNextPowerMethod())
        	{
        		for(Double alpha = params.getFirstAlpha(); alpha != null;
        		alpha = params.getNextAlpha())
        		{
        			for(Double sigmaScale = params.getFirstSigmaScale(); sigmaScale != null;
        			sigmaScale = params.getNextSigmaScale())
        			{
        				for(Double betaScale = params.getFirstBetaScale(); betaScale != null;
        				betaScale = params.getNextBetaScale())
        				{
        					for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
        					sampleSize = params.getNextSampleSize())
        					{                
        						Double quantile = params.getFirstQuantile(); 
        						do
        						{
        							// we only continue through this loop for quantile power
        							double power = getPowerByType(params);
        							if (method == PowerMethod.QUANTILE_POWER)
        								results.add(new GLMMPower(test, alpha, 
        									power, power, 
        									params.getDesignEssence().getTotalSampleSize(), 
        									betaScale, sigmaScale, method, quantile));
        							else
        								results.add(new GLMMPower(test, alpha, 
            									power, power, 
            									params.getDesignEssence().getTotalSampleSize(), 
            									betaScale, sigmaScale, method));
        							quantile = params.getNextQuantile();
        						}
        						while (method == PowerMethod.QUANTILE_POWER &&
        								quantile != null);
        					}
        				}
        			}
        		}
        	}
        }
        return results;
	}

	/**
	 * Find the best possible sample size to achieve a specified power or
	 * list of powers.  Sample size is determined with a bisection search
	 * 
	 * @see GLMMPowerParameters
	 * @see GLMMPower
	 * @param sampleSizeParams inputs to the sample size calculation
	 * @return list of calculated power values.
	 */
	@Override
	public List<Power> getSampleSize(PowerParameters sampleSizeParams)
	{
        GLMMPowerParameters params = (GLMMPowerParameters) sampleSizeParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);
        
        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();
        
        // calculate the power for either one or two tails
        for(GLMMPowerParameters.Test test = params.getFirstTest(); test != null;
        test = params.getNextTest())
        {           
        	for(GLMMPowerParameters.PowerMethod method = params.getFirstPowerMethod(); method != null;
        	method = params.getNextPowerMethod())
        	{
        		for(Double alpha = params.getFirstAlpha(); alpha != null;
        		alpha = params.getNextAlpha())
        		{
        			for(Double sigmaScale = params.getFirstSigmaScale(); sigmaScale != null;
        			sigmaScale = params.getNextSigmaScale())
        			{
        				for(Double betaScale = params.getFirstBetaScale(); betaScale != null;
        				betaScale = params.getNextBetaScale())
        				{
        					// we can't calculate a sample size for no difference between groups, so
        					// we ignore this case for now.
        					if (betaScale == 0) continue;  
        					for(Double power = params.getFirstPower(); power != null; 
        					power = params.getNextPower())
        					{        						
        						Double quantile = params.getFirstQuantile(); 
        						do
        						{
        							// we only continue through this loop for quantile power
        							try
        							{
        								SampleSize sampleSize = getSampleSize(params);
            							if (method == PowerMethod.QUANTILE_POWER)
            								results.add(new GLMMPower(test, alpha.doubleValue(), 
            										power.doubleValue(), sampleSize.actualPower, sampleSize.sampleSize, 
            										betaScale.doubleValue(), sigmaScale.doubleValue(), method, quantile));
            							else
            								results.add(new GLMMPower(test, alpha.doubleValue(), 
            										power.doubleValue(), sampleSize.actualPower, sampleSize.sampleSize, 
            										betaScale.doubleValue(), sigmaScale.doubleValue(), method));	
        							}
        							catch (Exception e)
        							{
        								System.err.println("Sample size failed: " + e.getMessage());
        								// TODO: 
        							}    							
        							quantile = params.getNextQuantile();      	
        						}
        						while (method == PowerMethod.QUANTILE_POWER &&
        								quantile != null);
        					}
        				}
        			}
        		}
        	}
        }

		return results;
	}

	/**
	 * Runs a simulation to determine power values for the given
	 * parameters.  
	 * 
	 * @see GLMMPowerParameters
	 * @see GLMMPower
	 * @param powerParams inputs to the simulation
	 * @param iterations number of iterations to perform in the simulation
	 * @return list of calculated power values.
	 */
	@Override
	public List<Power> getSimulatedPower(PowerParameters powerParams, int iterations)
	{
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);
        
        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();
        
        // simulate power for all variations of the study design
        for(GLMMPowerParameters.Test test = params.getFirstTest(); test != null;
        test = params.getNextTest())
        {
        	for(GLMMPowerParameters.PowerMethod method = params.getFirstPowerMethod(); method != null;
        	method = params.getNextPowerMethod())
        	{
        		for(Double alpha = params.getFirstAlpha(); alpha != null;
        		alpha = params.getNextAlpha())
        		{
        			for(Double sigmaScale = params.getFirstSigmaScale(); sigmaScale != null;
        			sigmaScale = params.getNextSigmaScale())
        			{
        				for(Double betaScale = params.getFirstBetaScale(); betaScale != null;
        				betaScale = params.getNextBetaScale())
        				{
        					for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
        					sampleSize = params.getNextSampleSize())
        					{
           						Double quantile = params.getFirstQuantile(); 
           						do
           						{
           							// simulate the power
           							double power = simulatePower(params, iterations);
           							// store the power result
        							if (method == PowerMethod.QUANTILE_POWER)
        								results.add(new GLMMPower(test, alpha, 
        									power, power, 
        									params.getDesignEssence().getTotalSampleSize(), 
        									betaScale, sigmaScale, method, quantile));
        							else
        								results.add(new GLMMPower(test, alpha, 
            									power, power, 
            									params.getDesignEssence().getTotalSampleSize(), 
            									betaScale, sigmaScale, method));
        							
        							quantile = params.getNextQuantile();      	
           						}
           						while (method == PowerMethod.QUANTILE_POWER &&
           								quantile != null);
        					}
        				}
        			}
        		}
        	}
        }
        
        return results;
	}

	/**
	 * Find the best possible effect size (i.e. scale factor for the beta matrix)  to achieve 
	 * a specified power or list of powers.  Effect size is determined with a bisection search
	 * 
	 * @see GLMMPowerParameters
	 * @see GLMMPower
	 * @param powerParams inputs to the effect size calculation
	 * @return list of calculated power values.
	 */
	@Override
	public List<Power> getDetectableDifference(PowerParameters powerParams)
	{
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);
        
        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();
        
        // calculate the power for all variations of the study design
        for(GLMMPowerParameters.Test test = params.getFirstTest(); test != null;
        test = params.getNextTest())
        {
        	for(GLMMPowerParameters.PowerMethod method = params.getFirstPowerMethod(); method != null;
        	method = params.getNextPowerMethod())
        	{
        		for(Double alpha = params.getFirstAlpha(); alpha != null;
        		alpha = params.getNextAlpha())
        		{
        			for(Double sigmaScale = params.getFirstSigmaScale(); sigmaScale != null;
        			sigmaScale = params.getNextSigmaScale())
        			{
        				for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null;
        				sampleSize = params.getNextSampleSize())
        				{
        					for(Double power = params.getFirstPower(); power != null; 
        					power = params.getNextPower())
        					{
        						Double quantile = params.getFirstQuantile(); 
        						do
        						{
        							try
        							{
        								DetectableDifference betaScale = getDetectableDifference(params);
        								results.add(new GLMMPower(test, alpha.doubleValue(), 
        										power.doubleValue(), betaScale.actualPower, sampleSize.intValue(), 
        										betaScale.betaScale, sigmaScale.doubleValue(), method));

               							if (method == PowerMethod.QUANTILE_POWER)
            								results.add(new GLMMPower(test, alpha.doubleValue(), 
            										power.doubleValue(), betaScale.actualPower, sampleSize.intValue(), 
            										betaScale.betaScale, sigmaScale.doubleValue(), method, quantile));
            							else
            								results.add(new GLMMPower(test, alpha.doubleValue(), 
            										power.doubleValue(), betaScale.actualPower, sampleSize.intValue(), 
            										betaScale.betaScale, sigmaScale.doubleValue(), method));
               							
            							quantile = params.getNextQuantile();      								
        							}
        							catch (Exception e)
        							{
        								System.err.println("Sample size failed: " + e.getMessage());
        								// TODO: 
        							}
        						}
        						while (method == PowerMethod.QUANTILE_POWER &&
        								quantile != null);
        					}
        				}
        			}
        		}
        	}
        }
        return results;
	}
	
	/**
	 * Perform any preliminary calculations / updates on the input matrices
	 * @param params
	 */
    private void initialize(GLMMPowerParameters params)
    {           	
        // update the sigma error if we have a baseline covariate
        DesignEssenceMatrix XEssence = params.getDesignEssence();
        int numRandom = 
        	(XEssence.getRandomMatrix() != null ? XEssence.getRandomMatrix().getColumnDimension() : 0);
        if (numRandom == 1)
        {
            RealMatrix sigmaG = params.getSigmaGaussianRandom();
            RealMatrix sigmaY = params.getSigmaOutcome();
            RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
            
            // set the sigma error matrix to [sigmaY - sigmaYG * sigmaG-1 * sigmaGY] 
            RealMatrix sigmaGY = sigmaYG.transpose();
            RealMatrix sigmaGInverse = new LUDecompositionImpl(sigmaG).getSolver().getInverse();
            params.setSigmaError(sigmaY.subtract(sigmaYG.multiply(sigmaGInverse.multiply(sigmaGY))));
            
            // calculate the betaG matrix and fill in the placeholder row for the random predictor
            FixedRandomMatrix beta = params.getBeta();
            beta.updateRandomMatrix(sigmaGInverse.multiply(sigmaGY));
        }
    }
	
    /**
     * Ensure that all required matrices are specified, and that conformance is correct
     * @param params
     * @throws IllegalArgumentException
     */
	protected void validateMatrices(GLMMPowerParameters params) throws IllegalArgumentException
	{
	       // convenience variables
        RealMatrix beta = params.getBeta().getCombinedMatrix();
        RealMatrix theta0 = params.getTheta();
        RealMatrix X = params.getDesign();
        DesignEssenceMatrix XEssence = params.getDesignEssence();
        int numRandom = 
        	(XEssence.getRandomMatrix() != null ? XEssence.getRandomMatrix().getColumnDimension() : 0);
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        
        // only allow at most 1 random predictor
        // TODO: handle multiple random predictors
        if (numRandom > 1)
            throw new IllegalArgumentException("Two many random predictors - at most 1 is allowed");

        // make sure all required matrices have been specified
        // note, we don't check U (within subject contrast), since it may be null in univariate cases
        if (beta == null) 
            throw new IllegalArgumentException("No beta (regression coefficients) matrix specified");
        if (X == null)
            throw new IllegalArgumentException("No design matrix specified");
        if (C == null)
            throw new IllegalArgumentException("No between subject contrast (C) matrix specified");
        if (theta0 == null)
            throw new IllegalArgumentException("No theta_null (null hypothesis) matrix specified");
        // create a default U if not specified
        if (U == null)
        {
            U = MatrixUtils.createRealIdentityMatrix(beta.getColumnDimension());
            params.setWithinSubjectContrast(U);
        }
        
        // different variance/covariance matrices are specified depending on the presence
        // of random covariate
        if (numRandom == 0)
        {
            if (sigmaE == null)
                throw new IllegalArgumentException("No sigma (error) matrix specified");
            if (!sigmaE.isSquare())
                throw new IllegalArgumentException("Sigma error matrix must be square");
            if (U.getRowDimension() != sigmaE.getRowDimension())
                throw new IllegalArgumentException("Within subject contrast does not conform with sigma matrix");
        }
        else if (numRandom == 1)
        {
            // TODO: FIX THIS make sure the test statistic is either HLT or UNIREP if there is a random
            // covariate (results not published for Wilk's Lambda or Pillai-Bartlett 
        	for(Test test = params.getFirstTest(); test != null; test = params.getNextTest())
        	{
                if (test != Test.HOTELLING_LAWLEY_TRACE && test != Test.UNIREP && test != Test.UNIREP_BOX &&
                		test != Test.UNIREP_GEISSER_GREENHOUSE && test != Test.UNIREP_HUYNH_FELDT)
                    throw new IllegalArgumentException("With a random covariate, only Hotelling-Lawley and Unirep test statistics are supported");
        	}

            if (sigmaG == null)
                throw new IllegalArgumentException("No variance/covariance matrix specified for gaussian predictors");
            if (sigmaY == null)
                throw new IllegalArgumentException("No variance/covariance matrix specified for response variables");
            if (sigmaYG == null)
                throw new IllegalArgumentException("No outcome / gaussian predictor covariance matrix specified");
            
            // check conformance
            if (U.getRowDimension() != sigmaY.getRowDimension())
                throw new IllegalArgumentException("Within subject contrast does not conform with sigma matrix");
            if (sigmaG.getRowDimension() != sigmaYG.getColumnDimension())
                throw new IllegalArgumentException("Outcome / Gaussian predictor covariance matrix does not conform with variance matrix for the gaussian predictor");
            if (!sigmaY.isSquare())
                throw new IllegalArgumentException("Variance/covariance matrix for response variables must be square");
            if (!sigmaG.isSquare())
                throw new IllegalArgumentException("Variance/covariance matrix for gaussian predictors must be square");
        }
        
        // check dimensionality 
        if (C.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Between subject contrast does not conform with beta matrix");
        if (beta.getColumnDimension() != U.getRowDimension())
            throw new IllegalArgumentException("Within subject contrast does not conform with beta matrix");
        if (X.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Design matrix does not conform with beta matrix");
        if (C.getRowDimension() > C.getColumnDimension())
            throw new IllegalArgumentException("Number of rows in between subject contrast must be less than or equal to the number of columns");
        if (U.getColumnDimension() > U.getRowDimension())
            throw new IllegalArgumentException("Number of columns in within subject contrast must be less than or equal to the number of rows");
        if (theta0.getRowDimension() != C.getRowDimension())
            throw new IllegalArgumentException("Number of rows in theta null must equal number of rows in between subject contrast");

        // check rank of the design matrix
        RealMatrix XFixed = XEssence.getFixedMatrix();
        int rankX = new SingularValueDecompositionImpl(XFixed).getRank();
        if (rankX != Math.min(XFixed.getColumnDimension(), XFixed.getRowDimension()))
            throw new IllegalArgumentException("Design matrix is not full rank");            

        // make sure design matrix is symmetric and positive definite
        // TODO: how to check this?		
	}
	
	/**
	 * Convenience function to determine which power method to use 
	 *
	 * @param params
	 * @return
	 */
	private double getPowerByType(GLMMPowerParameters params)
	{
        // calculate the power
        double power = Double.NaN;
        switch (params.getCurrentPowerMethod())
        {
        case QUANTILE_POWER:
            power = getQuantilePower(params);
            break;
        case UNCONDITIONAL_POWER:
            power = getUnconditionalPower(params);
            break;
        case CONDITIONAL_POWER:
        default:
            power = getConditionalPower(params);
            break;
        }
        
        return power;
	}
	
	/**
	 * Compute conditional power.  Conditional power is conditional on
	 * a single instance of the design matrix, and is most appropriate for
	 * designs with only categorical  predictors
	 * @param params
	 * @return
	 */
    private double getConditionalPower(GLMMPowerParameters params)
    {
    	// create a test
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic 
        // under the null hypothesis
        double nonCentralityParam = glmmTest.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE);

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));
    }
    
    /**
     * Calculate power by integrating over all possible values of the
     * non-centrality parameter.  Best used for designs with a 
     * baseline covariate
     * 
     * @param params
     * @return
     * @throws IllegalArgumentException
     */
    private double getUnconditionalPower(GLMMPowerParameters params)
    throws IllegalArgumentException
    {  		
    	// create a test
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
        // get the distribution of the noncentrality parameter
        NonCentralityDistribution nonCentralityDist = new NonCentralityDistribution(params, false);
        double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double h1 = nonCentralityDist.getH1();
        double h0 = nonCentralityDist.getH0();
        // integrate over all value of non-centrality parameter from h0 to h1
        UnconditionalPowerIntegrand integrand = 
            new UnconditionalPowerIntegrand(nonCentralityDist, Fcrit, ndf, ddf);
        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        try
        {
            // create a noncentral F dist with non-centrality of H1
            NonCentralFDistribution fdist = new NonCentralFDistribution(ndf, ddf, h1);
            double integralResult = 0;
            if (h1 != 0) integralResult = integrator.integrate(integrand, h0, h1);
            
            return (1 - fdist.cdf(Fcrit) - 0.5*integralResult);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to integrate over non-centrality parameter: " + e.getMessage());
        }        
    }
    
    /**
     * Calculate quantile power by determining a specified quantile
     * of the non-centrality distribution.
     * 
     * @param params
     * @return
     */
    private double getQuantilePower(GLMMPowerParameters params)
    {
    	// create a test
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        NonCentralityDistribution nonCentralityDist = 
            new NonCentralityDistribution(params, params.isNonCentralityCDFExact());
        double nonCentralityParam = nonCentralityDist.inverseCDF(params.getCurrentQuantile());
        
        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NonCentralFDistribution nonCentralFDist = new NonCentralFDistribution(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }

    /**
     * Calculate power for the general linear multivariate model based on
     * the input matrices.
     * 
     * For the mutivariate case, the following matrices must be specified
     * <ul>
     * <li>Essence design matrix - essence matrix for the model design
     * <li>Beta matrix - estimated regression coefficients matrix
     * <li>Sigma matrix - estimated errors matrix
     * <li>Theta0 matrix - estimated null hypothesis matrix 
     * <li>Between subject contrast matrix - defines comparisons between subjects
     * <li>Within subject contrast matrix - defines comparisons within subjects
     * (may be left null for univariate special case)
     * </ul>
     * 
     * @param params Container for input matrices
     * @return power
     */
    private SampleSize getSampleSize(GLMMPowerParameters params)
            throws IllegalArgumentException, MathException
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(params);
        
        // find the per group sample size 
        DesignEssenceMatrix essence = params.getDesignEssence();
        int upperBound = getSampleSizeUpperBound(params);
        int minSampleSize = essence.getMinimumSampleSize();
        int perGroupSampleSize = (int) Math.ceil(solver.solve(sampleSizeFunc, minSampleSize, upperBound));
        if (perGroupSampleSize < 0) 
            throw new MathException("Failed to calculate sample size");
        
        // calculate the total N and return the values
        params.setGroupSampleSize(perGroupSampleSize);
        return new SampleSize(essence.getTotalSampleSize(), getPowerByType(params));
    }

    /**
     * Determine the upper bound for the bisection search used in 
     * calculation of sample size
     * 
     * @param params
     * @return
     */
    private int getSampleSizeUpperBound(GLMMPowerParameters params)
    {
        double desiredPower = params.getCurrentPower();
        int upperBound = STARTING_SAMPLE_SIZE;
        
        for(double currentPower = 0.0; currentPower < desiredPower; upperBound *= 2)
        {
            params.setGroupSampleSize(upperBound); 
            currentPower = getPowerByType(params);
        }
        return upperBound;
    }

    /**
     * Perform a bisection search to determine effect size
     * @param params
     * @return
     * @throws IllegalArgumentException
     * @throws MathException
     */
    private DetectableDifference getDetectableDifference(GLMMPowerParameters params)
    throws IllegalArgumentException, MathException
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        DetectableDifferenceFunction ddFunc = new DetectableDifferenceFunction(params);

        // find the detectable difference (i.e. beta scale)
        int upperBound = (int) Math.ceil(getDetectableDifferenceUpperBound(params));
        double betaScale = solver.solve(ddFunc, 0, upperBound);
        if (betaScale < 0) 
            throw new MathException("Failed to calculate sample size");

        // calculate the total N and return the values
        params.scaleBeta(betaScale);
        return new DetectableDifference(betaScale, getPowerByType(params));
    }

    /**
     * Get the upper bound for the bisection search used to determine effect size
     * @param params
     * @return
     */
    private double getDetectableDifferenceUpperBound(GLMMPowerParameters params)
    {
        double desiredPower = params.getCurrentPower();
        int upperBound = STARTING_DETECTABLE_DIFFERENCE;

        for(double currentPower = 0.0; currentPower < desiredPower; upperBound *= 2)
        {
            params.scaleBeta((double) upperBound);
            currentPower = getPowerByType(params);
        }
        return upperBound;
    }
    
    /**
     * Simulate the error matrix in the Y = X * beta + e
     */
    private RealMatrix simulateError(Normal normalDist, RealMatrix error, int rows, int columns, RealMatrix sigma)
    throws IllegalArgumentException
    {        
        // build a matrix of random values from a standard normal
        // the number of rows = #subjects (rows) in the full design matrix
        // the number of columns = #outcome variables (i.e. columns in beta)
        for(int rowIndex = 0; rowIndex < rows; rowIndex++)
        {
            for(int columnIndex = 0; columnIndex < columns; columnIndex++)
            {
            	error.setEntry(rowIndex, columnIndex, normalDist.random()); 
            }
        }
        
        // take the square root of the sigma matrix via cholesky decomposition
        try
        {            
            RealMatrix sqrtMatrix = 
                new CholeskyDecompositionImpl(sigma, 
                        RELATIVE_SYMMETRY_THRESHOLD,
                        POSITIVITY_THRESHOLD).getLT();
            return error.multiply(sqrtMatrix); 
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
    }
    
    /**
     * Simulate power for the general linear multivariate model based on
     * the input matrices.
     * 
     * For the mutivariate case, the following matrices must be specified
     * <ul>
     * <li>Design matrix - essence matrix for the model design
     * <li>Beta matrix - estimated regression coefficients matrix
     * <li>Sigma matrix - estimated errors matrix
     * <li>Theta0 matrix - estimated null hypothesis matrix 
     * <li>Between subject contrast matrix - defines comparisons between subjects
     * <li>Within subject contrast matrix - defines comparisons within subjects
     * (may be left null for univariate special case)
     * </ul>
     * 
     * @param params Container for input matrices
     * @param iterations number of simulation samples/iterations
     * @return simulated power value
     */
    public double simulatePower(GLMMPowerParameters params, int iterations)
            throws IllegalArgumentException
    {        
    	// get total observations, N, and rank of design matrix
    	RealMatrix X = params.getDesign();
    	double N = X.getRowDimension();
    	double rankX = params.getDesignRank();

    	// create a normal distribution for generating random errors
    	Normal normalDist = new Normal();       
    	normalDist.setSeed(1234);
		int rejectionCount = 0;
		// create an error matrix here, so we don't have to reallocate every time
        Array2DRowRealMatrix randomNormals = 
        	new Array2DRowRealMatrix((int) N, params.getScaledBeta().getColumnDimension());
        
		// TODO: separate simulation function for quantile/unconditional vs conditional
    	if (params.getDesignEssence().hasRandom() && 
    			params.getCurrentPowerMethod() != PowerMethod.CONDITIONAL_POWER)
    	{
    		for(int gInstance = 0; gInstance < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; gInstance++)
    		{
    			X = params.getDesign(true); // force a new realization of the design matrix (i.e. a new covariate column)
    			for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
    			{
    				if (simulateAndFitModel(params, normalDist, randomNormals, N, rankX)) rejectionCount++;
    			}
    		}
            return ((double) rejectionCount) / (((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL) *
            		((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));

    	}
    	else 
    	{
    		// run the simulations
    		for(int i = 0; i < iterations; i++)
    		{
        		if (simulateAndFitModel(params, normalDist, randomNormals, N, rankX)) rejectionCount++;
    		}
            return ((double) rejectionCount) / ((double) iterations);
    	}
    }
    
    /**
     * Simulate the error matrix to generate a single realization of the data, then
     * fit the model and determine if the null hypothesis is rejected
     * 
     * @param params GLMM parameter set
     * @param normalDist normal distribution object for generating random errors
     * @param N total number of observations
     * @param rankX rank of the design matrix
     * @return
     */
    private boolean simulateAndFitModel(GLMMPowerParameters params, Normal normalDist, 
    		RealMatrix randomNormals, double N, double rankX)
    {
    	RealMatrix X = params.getDesign();
		RealMatrix scaledBeta = params.getScaledBeta();
		RealMatrix scaledSigma = params.getScaledSigmaError();

		RealMatrix error = simulateError(normalDist, randomNormals, X.getRowDimension(),
					scaledBeta.getColumnDimension(), scaledSigma);         

		// calculate simulated Y based on Y = X beta + error
		RealMatrix Ysim = (X.multiply(scaledBeta)).add(error);
		// calculate beta-Hat
		RealMatrix XtXInverse = params.getXtXInverse();
		RealMatrix betaHat = XtXInverse.multiply(X.transpose()).multiply(Ysim);
		// calculate Y-hat
		RealMatrix YHat = (X.multiply(betaHat));
		// calculate sigma-Hat
		RealMatrix Ydiff = Ysim.subtract(YHat);
		RealMatrix sigmaHat = (Ydiff.transpose().multiply(Ydiff)).scalarMultiply(((double) 1/(N - rankX)));    

		// set the scaled sigma / beta to the new values
		params.setScaledBeta(betaHat);
		params.setScaledSigmaError(sigmaHat);

		// calculate the observed F for the simulation
		GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
		double fobs = glmmTest.getObservedF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

		// get the p-value from a central F distribution
		double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
		double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

		FishersF fdist = new FishersF(ndf, ddf);
		double pvalue = 1 - fdist.cdf(fobs);

		// reset the parameter values to the original scaled beta/sigma
		params.setScaledBeta(scaledBeta);
		params.setScaledSigmaError(scaledSigma);
		
		// check if we reject the null hypothesis
		return (pvalue <= params.getCurrentAlpha());
    }
    
}
