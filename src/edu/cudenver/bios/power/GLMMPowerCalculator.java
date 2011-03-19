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

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.LUDecompositionImpl;
//import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;
import org.apache.commons.math.stat.StatUtils;

import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.MatrixUtils;
import edu.cudenver.bios.matrix.RandomErrorMatrix;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
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
    private static final int STARTING_SAMPLE_SIZE = 1000;
    private static final int STARTING_BETA_SCALE = 1;
    private static final int SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL = 1000;
    
    /**
     * container class for simulation info
     */
    private class SimulationFit
    {
    	public RealMatrix Y;
    	public RealMatrix Ypred;
    	public RealMatrix sigma;
    	public RealMatrix beta;
    	public double numeratorDF;
    	public double denominatorDF;
    	public double Fvalue;
    	public double Pvalue;
    	
    	public SimulationFit(double Pvalue, double Fvalue, 
    			double numeratorDF, double denominatorDF,
    			RealMatrix Y, RealMatrix Ypred, RealMatrix sigma, RealMatrix beta)
    	{
        	this.Y = Y;
        	this.sigma = sigma;
        	this.beta = beta;
        	this.Pvalue = Pvalue;
        	this.Fvalue = Fvalue;
        	this.numeratorDF = numeratorDF;
        	this.denominatorDF = denominatorDF;
    	}
    }
    
    /**
     * Function used with Apache's bisection solver to determine the 
     * per-group sample size which most closely achieves the desired power
     */
    private class SampleSizeFunction implements UnivariateRealFunction
    {
    	private GLMMTest glmmTest;
    	private NonCentralityDistribution nonCentralityDist;
    	private PowerMethod method;
    	private double alpha;
    	private double quantile;
    	private double targetPower;
    	
        public SampleSizeFunction(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
        		PowerMethod method, double targetPower, double alpha, double quantile)
        {
            this.glmmTest = glmmTest;
            this.nonCentralityDist = nonCentralityDist;
            this.method = method;
            this.targetPower = targetPower;
            this.alpha = alpha;
            this.quantile = quantile;            
        }
        
        public double value(double n)
        {
            try
            {
                glmmTest.setPerGroupSampleSize((int) n); 
                double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
                return targetPower - calculatedPower;
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
     * Function used with Apache's bisection solver to determine the 
     * per-group sample size which most closely achieves the desired power
     */
    private class DetectableDifferenceFunction implements UnivariateRealFunction
    {
    	private GLMMTest glmmTest;
    	private NonCentralityDistribution nonCentralityDist;
    	private PowerMethod method;
    	private double alpha;
    	private double quantile;
    	private double targetPower;
    	private FixedRandomMatrix beta;
        
        public DetectableDifferenceFunction(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
        		PowerMethod method, double targetPower, double alpha, double quantile, FixedRandomMatrix beta)
        {
        	this.glmmTest = glmmTest;
        	this.nonCentralityDist = nonCentralityDist;
        	this.method = method;
        	this.alpha = alpha;
        	this.quantile = quantile;
        	this.targetPower = targetPower;
        	this.beta = beta;
        }
        
        public double value(double betaScale)
        {
            try
            {
                glmmTest.setBeta(beta.scalarMultiply(betaScale, true));
                double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
                return targetPower - calculatedPower;
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
		for(GLMMTestFactory.Test test : params.getTestList())
		{            
			for(GLMMPowerParameters.PowerMethod method : params.getPowerMethodList())
			{
				for(Double alpha : params.getAlphaList()) 
				{
					for(Double sigmaScale : params.getSigmaScaleList())
					{
						for(Double betaScale : params.getBetaScaleList())
						{
							int qIdx = 0;
							do
							{
								double quantile = (params.getQuantileList().size() > 0 ?
										params.getQuantileList().get(qIdx) : Double.NaN);
								for(Integer sampleSize : params.getSampleSizeList())
								{           
									try
									{
										results.add(getPowerValue(params, test, method, alpha, sigmaScale, betaScale,
												sampleSize, quantile));
									}
									catch (Exception e)
									{
										System.err.println("Power calculation failed: " + e.getMessage());
									}
								}
								qIdx++;
							} while (method == PowerMethod.QUANTILE_POWER &&
									qIdx < params.getQuantileList().size());
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
        // calculate the power for all variations of the study design
        for(GLMMTestFactory.Test test : params.getTestList())
        {            
        	for(GLMMPowerParameters.PowerMethod method : params.getPowerMethodList())
        	{
        		for(Double alpha : params.getAlphaList()) 
        		{
        			for(Double sigmaScale : params.getSigmaScaleList())
        			{
        				for(Double betaScale : params.getBetaScaleList())
        				{
        					// we can't calculate a sample size for no difference between groups, so
        					// we ignore this case for now.
        					if (betaScale == 0) continue;  
        					int qIdx = 0;
        					do
        					{
        						double quantile = (params.getQuantileList().size() > 0 ?
        								params.getQuantileList().get(qIdx) : Double.NaN);
        						for(Double power : params.getPowerList())
        						{    
        							try
        							{       									
        								results.add(getSampleSizeValue(params, test, method, alpha, 
        										sigmaScale, betaScale, power, quantile));
        							}
        							catch (Exception e)
        							{
        								System.err.println("Sample size failed: " + e.getMessage());
        								// TODO: how to handle this?
        							}    							
        							qIdx++;
        						} 
        					} while (method == PowerMethod.QUANTILE_POWER &&
        								qIdx < params.getQuantileList().size());
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
	 * Note: quantile / unconditional power currently hard-coded to 1000
	 * random instances of the baseline covariate x1000 error simulations.
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
		// calculate the power for all variations of the study design
		for(GLMMTestFactory.Test test : params.getTestList())
		{            
			for(GLMMPowerParameters.PowerMethod method : params.getPowerMethodList())
			{
				for(Double alpha : params.getAlphaList()) 
				{
					for(Double sigmaScale : params.getSigmaScaleList())
					{
						for(Double betaScale : params.getBetaScaleList())
						{
							int qIdx = 0;
							do
							{
								double quantile = (params.getQuantileList().size() > 0 ?
										params.getQuantileList().get(qIdx) : Double.NaN);
								for(Integer sampleSize : params.getSampleSizeList())
								{       
           							// simulate the power
									try
									{
										results.add(getSimulatedPowerValue(params, test, method, alpha, 
												sigmaScale, betaScale, sampleSize, quantile, iterations));
									}
									catch (Exception e)
									{
										System.err.println();
									}
           						}
								qIdx++;
							} while (method == PowerMethod.QUANTILE_POWER &&
									qIdx < params.getQuantileList().size());
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
		for(GLMMTestFactory.Test test : params.getTestList())
		{            
			for(GLMMPowerParameters.PowerMethod method : params.getPowerMethodList())
			{
				for(Double alpha : params.getAlphaList()) 
				{
					for(Double sigmaScale : params.getSigmaScaleList())
					{
						for(Integer sampleSize : params.getSampleSizeList())
						{
							int qIdx = 0;
							do
							{
								double quantile = (params.getQuantileList().size() > 0 ?
										params.getQuantileList().get(qIdx) : Double.NaN);
								for(Double power : params.getPowerList())
								{    
        							try
        							{
        								results.add(getDetectableDifference(params,
        										test, method, alpha, sigmaScale, power, sampleSize, quantile));			
        							}
        							catch (Exception e)
        							{
        								System.err.println("Detectable difference failed: " + e.getMessage());
        								// TODO: 
        							}
        						}
								qIdx++;
							} while (method == PowerMethod.QUANTILE_POWER &&
									qIdx < params.getQuantileList().size());
        				}
        			}
        		}
        	}
        }
        return results;
	}
	
	/**
	 * Perform any preliminary calculations / updates on the input matrices
	 * @param params input parameters
	 */
    private void initialize(GLMMPowerParameters params)
    {           	
        // update the sigma error if we have a baseline covariate
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        int numRandom = 
        	(sigmaG != null ? sigmaG.getRowDimension() : 0);
        if (numRandom == 1)
        {           
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
     * @param params GLMM input parameters
     * @throws IllegalArgumentException
     */
	protected void validateMatrices(GLMMPowerParameters params) throws IllegalArgumentException
	{
	       // convenience variables
        RealMatrix beta = params.getBeta().getCombinedMatrix();
        RealMatrix theta0 = params.getTheta();
        RealMatrix XEssence = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        int numRandom = 
        	(sigmaG != null ? sigmaG.getRowDimension() : 0);
        
        // only allow at most 1 random predictor
        // TODO: handle multiple random predictors
        if (numRandom > 1)
            throw new IllegalArgumentException("Two many random predictors - at most 1 is allowed");

        // make sure all required matrices have been specified
        // note, we don't check U (within subject contrast), since it may be null in univariate cases
        if (beta == null) 
            throw new IllegalArgumentException("No beta (regression coefficients) matrix specified");
        if (XEssence == null)
            throw new IllegalArgumentException("No design essence matrix specified");
        if (C == null)
            throw new IllegalArgumentException("No between subject contrast (C) matrix specified");
        if (theta0 == null)
            throw new IllegalArgumentException("No theta_null (null hypothesis) matrix specified");
        // create a default U if not specified
        if (U == null)
        {
            U = org.apache.commons.math.linear.MatrixUtils.createRealIdentityMatrix(beta.getColumnDimension());
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
            // covariate (results not published for Wilk's Lambda or Pillai-Bartlett 
        	for(Test test : params.getTestList())
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
        if (XEssence.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Design matrix does not conform with beta matrix");
        if (C.getRowDimension() > C.getColumnDimension())
            throw new IllegalArgumentException("Number of rows in between subject contrast must be less than or equal to the number of columns");
        if (U.getColumnDimension() > U.getRowDimension())
            throw new IllegalArgumentException("Number of columns in within subject contrast must be less than or equal to the number of rows");
        if (theta0.getRowDimension() != C.getRowDimension())
            throw new IllegalArgumentException("Number of rows in theta null must equal number of rows in between subject contrast");

        // check rank of the design matrix
        int rankX = new SingularValueDecompositionImpl(XEssence).getRank();
        if (rankX != Math.min(XEssence.getColumnDimension(), XEssence.getRowDimension()))
            throw new IllegalArgumentException("Design matrix is not full rank");            

        // make sure design matrix is symmetric and positive definite
        // TODO: how to check this?		
	}
		
	/**
	 * Compute conditional power.  Conditional power is conditional on
	 * a single instance of the design matrix, and is most appropriate for
	 * designs with only categorical  predictors
	 * @param params GLMN input parameters
	 * @return conditional power
	 */
    private double getConditionalPower(GLMMTest glmmTest, double alpha)
    {        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
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
     * @param params GLMM input parameters
     * @return unconditional power
     * @throws IllegalArgumentException
     */
    private double getUnconditionalPower(GLMMTest glmmTest, 
    		NonCentralityDistribution nonCentralityDist, double alpha)
    throws IllegalArgumentException
    {  		        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
        // get the distribution of the noncentrality parameter
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
     * @param params GLMM input parameters
     * @return quantile power
     */
    private double getQuantilePower(GLMMTest glmmTest, 
    		NonCentralityDistribution nonCentralityDist, double alpha, double quantile)
    {        
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        double nonCentralityParam = nonCentralityDist.inverseCDF(quantile);
        
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
     *  Find the sample size which achieves the desired power(s) 
     *  specified in the input parameters.  Uses a bisection search.
     * 
     * @param params GLMM input parameters
     * @return sample size
     */
    private GLMMPower getSampleSizeValue(GLMMPowerParameters params,
    		Test test, PowerMethod method, double alpha, 
    		double sigmaScale, double betaScale, double targetPower, double quantile)
            throws IllegalArgumentException, MathException
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

    	RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
    	RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
    	GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(test, 
				params.getFApproximationMethod(test),
				params.getUnivariateCdfMethod(test), 
				params.getDesignEssence(), 
				params.getXtXInverse(), 
				STARTING_SAMPLE_SIZE,
				params.getDesignRank(), 
				params.getBetweenSubjectContrast().getCombinedMatrix(), 
				params.getWithinSubjectContrast(), 
				params.getTheta(), 
				scaledBeta, scaledSigmaError);

    	NonCentralityDistribution nonCentralityDist = null;
    	if (method != PowerMethod.CONDITIONAL_POWER)
    	{
    		nonCentralityDist = new NonCentralityDistribution(test, 
    				params.getDesignEssence(), 
    				params.getXtXInverse(), 
    				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), STARTING_SAMPLE_SIZE),
    				params.getBetweenSubjectContrast(),
    				params.getWithinSubjectContrast(),
    				params.getTheta(),
    				scaledBeta, scaledSigmaError,
    				params.getSigmaGaussianRandom(),
    				params.isNonCentralityCDFExact());
    	}
        
    	// create a bisection search function to find the best per group sample size
        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(glmmTest, nonCentralityDist,
        		method, targetPower, alpha, quantile);
        
        // find the per group sample size 
        int upperBound = getSampleSizeUpperBound(glmmTest, nonCentralityDist,
        		method, targetPower, alpha, quantile);
        int perGroupSampleSize = (int) Math.ceil(solver.solve(sampleSizeFunc, 2, upperBound));
        if (perGroupSampleSize < 0) 
            throw new MathException("Failed to calculate sample size");
        // get the actual power associated with this per group sample size
        glmmTest.setPerGroupSampleSize(perGroupSampleSize);
        double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
		GLMMPowerConfidenceInterval ci = null;
		if (params.getConfidenceIntervalType() != 
			ConfidenceIntervalType.NONE)
		{        							
			ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
					params.getAlphaLowerConfidenceLimit(), 
					params.getAlphaUpperConfidenceLimit(), 
					params.getSampleSizeForEstimates(),
					params.getDesignMatrixRankForEstimates(), 
					alpha, glmmTest);
		}
        
        return new GLMMPower(test, alpha, targetPower, calculatedPower, 
        		MatrixUtils.getTotalSampleSize(params.getDesignEssence(), perGroupSampleSize), 
        		betaScale, sigmaScale, method, quantile, ci);
    }

    /**
     * Determine the upper bound for the bisection search used in 
     * calculation of sample size
     * 
     * @param params GLMM input parameters
     * @return upper bound on sample size to achieve desired power
     */
    private int getSampleSizeUpperBound(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
    		PowerMethod method, double targetPower, double alpha, double quantile)
    {
        int upperBound = STARTING_SAMPLE_SIZE;
        
        for(double currentPower = 0.0; currentPower < targetPower; upperBound *= 2)
        {
            glmmTest.setPerGroupSampleSize(upperBound); 
            currentPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
        }
        return upperBound;
    }

    /**
     * Perform a bisection search to determine effect size
     * @param params GLMM input parameters
     * @return detectable difference object
     * @throws IllegalArgumentException
     * @throws MathException
     */
    private GLMMPower getDetectableDifference(GLMMPowerParameters params,
    		Test test, PowerMethod method, double alpha, 
    		double sigmaScale, double targetPower, int sampleSize, double quantile)
    throws IllegalArgumentException, MathException
    {
    	RealMatrix scaledBeta = params.getBeta().scalarMultiply(STARTING_BETA_SCALE, true);
    	RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
    	GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(test, 
				params.getFApproximationMethod(test),
				params.getUnivariateCdfMethod(test), 
				params.getDesignEssence(), 
				params.getXtXInverse(), 
				sampleSize,
				params.getDesignRank(), 
				params.getBetweenSubjectContrast().getCombinedMatrix(), 
				params.getWithinSubjectContrast(), 
				params.getTheta(), 
				scaledBeta, scaledSigmaError);

    	NonCentralityDistribution nonCentralityDist = null;
    	if (method != PowerMethod.CONDITIONAL_POWER)
    	{
    		nonCentralityDist = new NonCentralityDistribution(test, 
    				params.getDesignEssence(), 
    				params.getXtXInverse(), 
    				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
    				params.getBetweenSubjectContrast(),
    				params.getWithinSubjectContrast(),
    				params.getTheta(),
    				scaledBeta, scaledSigmaError,
    				params.getSigmaGaussianRandom(),
    				params.isNonCentralityCDFExact());
    	}
    	
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();
        
        DetectableDifferenceFunction ddFunc = new DetectableDifferenceFunction(glmmTest, 
        		nonCentralityDist, method, targetPower, alpha, quantile, params.getBeta());

        // find the detectable difference (i.e. beta scale)
        int upperBound = (int) Math.ceil(getDetectableDifferenceUpperBound(glmmTest, 
        		nonCentralityDist, params.getBeta(),method, targetPower, alpha, quantile));
        double betaScale = solver.solve(ddFunc, 0, upperBound);
        if (betaScale < 0) 
            throw new MathException("Failed to calculate sample size");

        // calculate actual power associated with this beta scale
        glmmTest.setBeta(params.getBeta().scalarMultiply(betaScale, true));
        double calculatedPower = getPowerByType(glmmTest, nonCentralityDist, method,  alpha, quantile);
        
        // get a confidence interval if requested
		GLMMPowerConfidenceInterval ci = null;
		if (params.getConfidenceIntervalType() != 
			ConfidenceIntervalType.NONE)
		{        							
			ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
					params.getAlphaLowerConfidenceLimit(), 
					params.getAlphaUpperConfidenceLimit(), 
					params.getSampleSizeForEstimates(),
					params.getDesignMatrixRankForEstimates(), 
					alpha, glmmTest);
		}
        
		return new GLMMPower(test, alpha, targetPower, calculatedPower, 
				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize), betaScale, 
				sigmaScale, method, quantile, ci);
    }

    /**
     * Get the upper bound for the bisection search used to determine effect size
     * @param params GLMM input parameters
     * @return upper bound on beta scale
     */
    private double getDetectableDifferenceUpperBound(GLMMTest glmmTest, 
    		NonCentralityDistribution nonCentralityDist, FixedRandomMatrix beta,
    		PowerMethod method, double targetPower, double alpha, double quantile)
    {
        int upperBound = STARTING_BETA_SCALE;

        for(double currentPower = 0.0; currentPower < targetPower; upperBound *= 2)
        {
            glmmTest.setBeta(beta.scalarMultiply(upperBound, true)); 
            currentPower = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
        }
        return upperBound;
    }
    
    
//    /**
//     * Simulate power for the general linear multivariate model based on
//     * the input matrices.
//     * 
//     * @param params Container for input matrices
//     * @param iterations number of simulation samples/iterations
//     * @return simulated power value
//     */
//    public double simulatePower(GLMMPowerParameters params, int iterations)
//            throws IllegalArgumentException
//    {        
//    	// get total observations, N, and rank of design matrix
//    	RealMatrix XEssence = params.getDesignEssence();
//    	double N = XEssence.getRowDimension()*params.getCurrentSampleSize();
//    	double rankX = params.getDesignRank();
//
//    	// create a normal distribution for generating random errors
//    	Normal normalDist = new Normal();       
//    	normalDist.setSeed(1234);
//		int rejectionCount = 0;
//		// create an error matrix here, so we don't have to reallocate every time
//        Array2DRowRealMatrix randomNormals = 
//        	new Array2DRowRealMatrix((int) N, params.getScaledBeta().getColumnDimension());
//        
//    	if (params.getSigmaGaussianRandom() != null && 
//    			params.getCurrentPowerMethod() != PowerMethod.CONDITIONAL_POWER)
//    	{
//    		double[] powerValues = new double[SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL];
//    		
//    		for(int gInstance = 0; gInstance < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; gInstance++)
//    		{
//    			// force a new realization of the design matrix (i.e. a new covariate column)
//    			RealMatrix X = getFullDesignMatrix(params.getDesignEssence(), 
//    					params.getSigmaGaussianRandom(),	perGroupSize);
//    			rejectionCount = 0;
//    			for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
//    			{
//    				GLMMPowerParameters params, Test test,
//    	    		RealMatrix X, RealMatrix XtXinverse, 
//    	    		RealMatrix scaledBeta, 
//    				RandomErrorMatrix randomErrors, 
//    				int perGroupSampleSize, double alpha
//    				
//    				
//    				if (simulateAndFitModel(params, test, X, null, scaledBeta, random normalDist, randomNormals, N, rankX)) rejectionCount++;
//    			}
//    			powerValues[gInstance] = (((double) rejectionCount) / 
//    					((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
//    		}
//    		
//    		switch (params.getCurrentPowerMethod())
//    		{
//    		case UNCONDITIONAL_POWER:
//    			return StatUtils.mean(powerValues);
//    		case QUANTILE_POWER:
//    			return StatUtils.percentile(powerValues, params.getCurrentQuantile());
//    		default:
//    			throw new IllegalArgumentException("Unknown power method");
//    		}
//    	}
//    	else 
//    	{
//    		// run the simulations
//    		for(int i = 0; i < iterations; i++)
//    		{
//        		if (simulateAndFitModel(params, normalDist, randomNormals, N, rankX)) rejectionCount++;
//    		}
//            return ((double) rejectionCount) / ((double) iterations);
//    	}
//    }
    
    /**
     * Simulate the error matrix to generate a single realization of the data, then
     * fit the model and determine if the null hypothesis is rejected
     * 
     * @param params GLMM parameter set
     * @param normalDist normal distribution object for generating random errors
     * @param N total number of observations
     * @param rankX rank of the design matrix
     * @return true if null is rejected, false otherwise
     */
    private SimulationFit simulateAndFitModel(GLMMPowerParameters params, Test test,
    		RealMatrix X, RealMatrix XtXinverse, 
    		RealMatrix scaledBeta, 
			RandomErrorMatrix randomErrors, 
			int perGroupSampleSize, double alpha)
    {
		RealMatrix error = randomErrors.random();
		int N = X.getRowDimension();
		int rank = params.getDesignRank();
		
		// calculate simulated Y based on Y = X beta + error
		RealMatrix Ysim = (X.multiply(scaledBeta)).add(error);
		// calculate beta-Hat
		if (XtXinverse == null)
		{
			XtXinverse = params.getXtXInverse();
		}
		RealMatrix betaHat = XtXinverse.multiply(X.transpose()).multiply(Ysim);
		// calculate Y-hat
		RealMatrix YHat = (X.multiply(betaHat));
		// calculate sigma-Hat
		RealMatrix Ydiff = Ysim.subtract(YHat);
		RealMatrix sigmaHat = (Ydiff.transpose().multiply(Ydiff)).scalarMultiply(((double) 1/(N - rank)));    

		// build a test object for the simulated matrices
    	GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(test, 
				params.getFApproximationMethod(test),
				params.getUnivariateCdfMethod(test), 
				params.getDesignEssence(), 
				XtXinverse, 
				perGroupSampleSize,
				params.getDesignRank(), 
				params.getBetweenSubjectContrast().getCombinedMatrix(), 
				params.getWithinSubjectContrast(), 
				params.getTheta(), 
				betaHat, sigmaHat);
		// calculate the observed F for the simulation
		double fobs = glmmTest.getObservedF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

		// get the p-value from a central F distribution
		double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
		double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

		FishersF fdist = new FishersF(ndf, ddf);
		double pvalue = 1 - fdist.cdf(fobs);
		
		// return the fit information for the simulated matrices
		return new SimulationFit(pvalue, fobs, ndf, ddf, Ysim, YHat, sigmaHat, betaHat);
    }
    
    /**
     * Returns a power sample for each combination of parameters for a 
     * design with a random covariate (GLMM(F,g)).  Currently produces 
     * a sample of size 1000
     * 
     * @param params power parameters
     * @param size size of the sample
     * @return list of power samples
     * @throws IllegalArgumentException
     */
    public List<double[]> getSimulatedPowerSample(GLMMPowerParameters params, int size)
    throws IllegalArgumentException
    {		    
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);

    	if (params == null || params.getSigmaGaussianRandom() != null)
    		throw new IllegalArgumentException("Power samples can only be generated for designs with random predictors");
    	if (size <= 0) throw new IllegalArgumentException("Iterations must be positive");
    	
        // list of power results
        ArrayList<double[]> results = new ArrayList<double[]>();
        
		for(GLMMTestFactory.Test test : params.getTestList())
		{            
			for(GLMMPowerParameters.PowerMethod method : params.getPowerMethodList())
			{
				// only generate samples for quantile or unconditional power
				if (method == PowerMethod.CONDITIONAL_POWER) continue;
				for(Double alpha : params.getAlphaList()) 
				{
					for(Double sigmaScale : params.getSigmaScaleList())
					{
						for(Double betaScale : params.getBetaScaleList())
						{
							int qIdx = 0;
							do
							{
								double quantile = (params.getQuantileList().size() > 0 ?
										params.getQuantileList().get(qIdx) : Double.NaN);
								for(Integer sampleSize : params.getSampleSizeList())
								{       
           							results.add(getSimulatedPowerSampleValue(params,
           						    		test, method, alpha, sigmaScale, betaScale, sampleSize, 
           						    		quantile, size));
    								qIdx++;
								}
        					} while (method == PowerMethod.QUANTILE_POWER &&
									qIdx < params.getQuantileList().size());
        				}
        			}
        		}
        	}
        }
        
        return results;
    }
    
    /**
     * Generate a sample of powers for a design with a random covariate
     * @param params power parameters
     * @param iterations size of the sample
     * @return
     */
    private double[] getSimulatedPowerSampleValue(GLMMPowerParameters params,
    		Test test, PowerMethod method, double alpha, 
    		double sigmaScale, double betaScale, int sampleSize, double quantile, int iterations)
    {
    	// scale beta and sigma
    	RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
    	RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
    	// get random errors
    	RandomErrorMatrix randomErrors = 
    		new RandomErrorMatrix(MatrixUtils.getTotalSampleSize(params.getDesignEssence(), 
    				sampleSize), scaledBeta.getColumnDimension(), scaledSigmaError);
    	
		double[] powerValues = new double[iterations];

		for(int gInstance = 0; gInstance < iterations; gInstance++)
		{
			// force a new realization of the design matrix (i.e. a new covariate column)
    		RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(), 
    				params.getSigmaGaussianRandom(), sampleSize);
    		int rejectionCount = 0;
			for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
			{
				SimulationFit fit = 
					simulateAndFitModel(params, test, X, null, scaledBeta, randomErrors, 
						sampleSize, alpha);
				if (fit.Pvalue <= alpha) rejectionCount++;
			}
			powerValues[gInstance] = (((double) rejectionCount) / 
					((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
		}
    	
    	return powerValues;
    }

	/**
	 * Convenience function to determine which power method to use 
	 *
	 * @param params GLMM input parameters
	 * @return conditional/quantile/unconditional power
	 */
    private GLMMPower getPowerValue(GLMMPowerParameters params,
    		Test test, PowerMethod method, double alpha, 
    		double sigmaScale, double betaScale, int sampleSize, double quantile)
    {
    	RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
    	RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
    	GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(test, 
				params.getFApproximationMethod(test),
				params.getUnivariateCdfMethod(test), 
				params.getDesignEssence(), 
				params.getXtXInverse(), 
				sampleSize,
				params.getDesignRank(), 
				params.getBetweenSubjectContrast().getCombinedMatrix(), 
				params.getWithinSubjectContrast(), 
				params.getTheta(), 
				scaledBeta, scaledSigmaError);

    	NonCentralityDistribution nonCentralityDist = null;
    	if (method != PowerMethod.CONDITIONAL_POWER)
    	{
    		nonCentralityDist = new NonCentralityDistribution(test, 
    				params.getDesignEssence(), 
    				params.getXtXInverse(), 
    				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize),
    				params.getBetweenSubjectContrast(),
    				params.getWithinSubjectContrast(),
    				params.getTheta(),
    				scaledBeta, scaledSigmaError,
    				params.getSigmaGaussianRandom(),
    				params.isNonCentralityCDFExact());
    	}

        // calculate the power
        double power = getPowerByType(glmmTest, nonCentralityDist, method, alpha, quantile);
        switch (method)
        {
        case QUANTILE_POWER:
            power = getQuantilePower(glmmTest, nonCentralityDist, alpha, quantile);
            break;
        case UNCONDITIONAL_POWER:
            power = getUnconditionalPower(glmmTest, nonCentralityDist, alpha);
            break;
        case CONDITIONAL_POWER:
        default:
            power = getConditionalPower(glmmTest, alpha);
            break;
        }
    	
		GLMMPowerConfidenceInterval ci = null;
		if (params.getConfidenceIntervalType() != 
			ConfidenceIntervalType.NONE)
		{        							
			ci = new GLMMPowerConfidenceInterval(params.getConfidenceIntervalType(),
					params.getAlphaLowerConfidenceLimit(), 
					params.getAlphaUpperConfidenceLimit(), 
					params.getSampleSizeForEstimates(),
					params.getDesignMatrixRankForEstimates(), 
					alpha, glmmTest);
		}
		
		return new GLMMPower(test, alpha, power, power, 
				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize), betaScale, 
				sigmaScale, method, quantile, ci);
    }
    
    private GLMMPower getSimulatedPowerValue(GLMMPowerParameters params,
    		Test test, PowerMethod method, double alpha, 
    		double sigmaScale, double betaScale, int sampleSize, double quantile,
    		int iterations)
    {
    	// scale beta and sigma
    	RealMatrix scaledBeta = params.getBeta().scalarMultiply(betaScale, true);
    	RealMatrix scaledSigmaError = params.getSigmaError().scalarMultiply(sigmaScale);
    	// get random errors
    	RandomErrorMatrix randomErrors = 
    		new RandomErrorMatrix(MatrixUtils.getTotalSampleSize(params.getDesignEssence(), 
    				sampleSize), scaledBeta.getColumnDimension(), scaledSigmaError);
    	
    	double power = Double.NaN;
    	int rejectionCount = 0;
    	if (params.getSigmaGaussianRandom() != null && 
    			method != PowerMethod.CONDITIONAL_POWER)
    	{
    		double[] powerValues = new double[SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL];
    		
    		for(int gInstance = 0; gInstance < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; gInstance++)
    		{
    			// force a new realization of the design matrix (i.e. a new covariate column)
        		RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(), 
        				params.getSigmaGaussianRandom(), sampleSize);
    			for(int i = 0; i < SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL; i++)
    			{
    				SimulationFit fit = 
    					simulateAndFitModel(params, test, X, null, scaledBeta, randomErrors, 
    						sampleSize, alpha);
    				if (fit.Pvalue <= alpha) rejectionCount++;
    			}
    			powerValues[gInstance] = (((double) rejectionCount) / 
    					((double) SIMULATION_ITERATIONS_QUANTILE_UNCONDITIONAL));
    		}
    		
    		switch (method)
    		{
    		case UNCONDITIONAL_POWER:
    			power = StatUtils.mean(powerValues);
    		case QUANTILE_POWER:
    			power= StatUtils.percentile(powerValues, quantile);
    		default:
    			throw new IllegalArgumentException("Unknown power method");
    		}
    	}
    	else
    	{
    		// GLMM(F) design, conditional power
    		// run the simulations
    		RealMatrix X = MatrixUtils.getFullDesignMatrix(params.getDesignEssence(), sampleSize);
    		for(int i = 0; i < iterations; i++)
    		{
				SimulationFit fit = 
					simulateAndFitModel(params, test, X, params.getXtXInverse().scalarMultiply(1/(double)sampleSize), 
							scaledBeta, randomErrors, 
						sampleSize, alpha);
				if (fit.Pvalue <= alpha) rejectionCount++;
    		}
            power =  ((double) rejectionCount) / ((double) iterations);
    	}
	    	
		return new GLMMPower(test, alpha, power, power, 
				MatrixUtils.getTotalSampleSize(params.getDesignEssence(), sampleSize), betaScale, 
				sigmaScale, method, quantile, null);
    }
    
    /**
     * Get an individual power instance
     * @param glmmTest
     * @param nonCentralityDist
     * @param method
     * @param targetPower
     * @param alpha
     * @param quantile
     * @return
     */
    private double getPowerByType(GLMMTest glmmTest, NonCentralityDistribution nonCentralityDist,
    		PowerMethod method, double alpha, double quantile)
    {
        // calculate the power
        double power = Double.NaN;
        switch (method)
        {
        case QUANTILE_POWER:
            power = getQuantilePower(glmmTest, nonCentralityDist, alpha, quantile);
            break;
        case UNCONDITIONAL_POWER:
            power = getUnconditionalPower(glmmTest, nonCentralityDist, alpha);
            break;
        case CONDITIONAL_POWER:
        default:
            power = getConditionalPower(glmmTest, alpha);
            break;
        }
        return power;
    }
}
