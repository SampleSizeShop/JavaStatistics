package edu.cudenver.bios.power;

import java.security.InvalidAlgorithmParameterException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;

import jsc.distributions.StudentsT;

import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;

public class OneSampleStudentsTPowerCalculator implements PowerCalculator
{
    private static final int MAX_SAMPLE_SIZE = 100000;
    private static final int MIN_SAMPLE_SIZE =  2; // need df > 0
    
    // function to be used with apache's built-in bisection solver
    private class SampleSizeFunction implements UnivariateRealFunction
    {
        double mu0;
        double muA;
        double sigma;
        double alpha;
        double power;
        boolean oneTailed = false;
        
        public SampleSizeFunction(double mu0, double muA, double sigma,
                double alpha, double power, boolean oneTailed)
        {
            this.mu0 = mu0;
            this.muA = muA;
            this.sigma = sigma;
            this.alpha = alpha;
            this.power = power;
            this.oneTailed = oneTailed;
        }
        
        public double value(double n)
        {
            try
            {
                // create a t distribution with the specified degrees of freedom
                double df = n -1;
                StudentsT tdist = new StudentsT(df);

                if (oneTailed)
                {
                    double tAlpha = tdist.inverseCdf(1-alpha);
                    double tBeta = tdist.inverseCdf(power);
                    double root = ((sigma*sigma*(tBeta + tAlpha)*(tBeta + tAlpha))/
                            ((mu0 - muA)*(mu0-muA))) - n;
                    return root;
                }
                else
                {
                    double tAlpha = tdist.inverseCdf(1-alpha/2);
                    double tBeta = tdist.inverseCdf(power);
                    double root = ((sigma*sigma*(tBeta + tAlpha)*(tBeta + tAlpha))/
                            ((mu0 - muA)*(mu0-muA))) - n;
                    return root;
                }
            }
            catch (Exception e)
            {   
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a large number to prevent BisectionSolver from using
                // the n which caused to exception as a root
                return MAX_SAMPLE_SIZE;  
            }
        }
    }
    
	@Override
	public List<Power> getDetectableDifference(PowerParameters params)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Power> getPower(PowerParameters params)
	{
        OneSampleStudentsTPowerParameters studentsTParams = (OneSampleStudentsTPowerParameters) params;
        
        ArrayList<Power> results = new ArrayList<Power>();
        	
        // calculate the power for either one or two tails
        for(Double alpha: studentsTParams.getAlphaList())
        {
        	for(Double sigma: studentsTParams.getSigmaList())
        	{
        		for(OneSampleStudentsTPowerParameters.MeanPair means: studentsTParams.getMeansList())
        		{
        			for(Integer sampleSize: studentsTParams.getSampleSizeList())
        			{
        	        if (studentsTParams.isOneTailed())
        	        {
        	            results.add(calculateOneTailedPower(means.mu0, means.muA, 
        	            		alpha.doubleValue(), sigma.doubleValue(), sampleSize.intValue()));
        	        }
        	        else
        	        {
        	        	results.add(calculateTwoTailedPower(means.mu0, means.muA, 
        	            		alpha.doubleValue(), sigma.doubleValue(), sampleSize.intValue()));
        	        }
        			}
        		}
        	}
        }

        return results;
	}

	@Override
	public List<Power> getSampleSize(PowerParameters params)
	{
        OneSampleStudentsTPowerParameters studentsTParams = (OneSampleStudentsTPowerParameters) params;
        
        ArrayList<Power> results = new ArrayList<Power>();
        	
        // calculate the power for either one or two tails
        for(Double alpha: studentsTParams.getAlphaList())
        {
        	for(Double sigma: studentsTParams.getSigmaList())
        	{
        		for(OneSampleStudentsTPowerParameters.MeanPair means: studentsTParams.getMeansList())
        		{
        			for(Double power: studentsTParams.getPowerList())
        			{
        				if (studentsTParams.isOneTailed())
        				{
        					results.add(calculateOneTailedSampleSize(means.mu0, means.muA, 
        							alpha.doubleValue(), sigma.doubleValue(), power.doubleValue()));
        				}
        				else
        				{
        					results.add(calculateTwoTailedSampleSize(means.mu0, means.muA, 
        							alpha.doubleValue(), sigma.doubleValue(), power.doubleValue()));
        				}
        			}
        		}
        	}
        }

        return results;
	}

	@Override
	public List<Power> getSimulatedPower(PowerParameters params, int iterations)
	{
        OneSampleStudentsTPowerParameters studentsTParams = (OneSampleStudentsTPowerParameters) params;
        
        ArrayList<Power> results = new ArrayList<Power>();
        	
        // calculate the power for either one or two tails
        for(Double alpha: studentsTParams.getAlphaList())
        {
        	for(Double sigma: studentsTParams.getSigmaList())
        	{
        		for(OneSampleStudentsTPowerParameters.MeanPair means: studentsTParams.getMeansList())
        		{
        			for(Integer sampleSize: studentsTParams.getSampleSizeList())
        			{
        	        if (studentsTParams.isOneTailed())
        	        {
        	            results.add(simulateOneTailedPower(means.mu0, means.muA, 
        	            		alpha.doubleValue(), sigma.doubleValue(), sampleSize.intValue(), iterations));
        	        }
        	        else
        	        {
        	        	results.add(simulateTwoTailedPower(means.mu0, means.muA, 
        	            		alpha.doubleValue(), sigma.doubleValue(), sampleSize.intValue(), iterations));
        	        }
        			}
        		}
        	}
        }

        return results;
	}

	   /*
     * Calculate a one-tailed power
     */
    protected OneSampleStudentsTPower calculateOneTailedPower(double mu0, double muA, double alpha,
            double sigma, int sampleSize) throws IllegalArgumentException
    {
        /* calculate the degrees of freedom */
        int df = sampleSize - 1;
        /* get the critical t for specified alpha level in H0 distribution (central T) */
        StudentsT tdist = new StudentsT(df);
        double tAlpha = tdist.inverseCdf(alpha);
        
        double power = (tdist.cdf(tAlpha + Math.abs(mu0 - muA)*Math.sqrt(sampleSize)/sigma));
        
        return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
    }

    /*
     * Calculate a two-tailed power
     */
    protected OneSampleStudentsTPower calculateTwoTailedPower(double mu0, double muA, double alpha,
            double sigma, int sampleSize) throws IllegalArgumentException
    {
        /* calculate the degrees of freedom */
        int df = sampleSize - 1;
        double meanDiff = Math.abs(mu0 - muA);
        /* get the critical t's for specified alpha level in H0 distribution (central T)*/
        StudentsT tdist = new StudentsT(df);
        double tAlphaLower = tdist.inverseCdf(alpha/2);
        double tAlphaUpper = tdist.inverseCdf(1-alpha/2);
        
        /* return the calculated power */
        double power = (tdist.cdf(tAlphaLower + meanDiff*Math.sqrt(sampleSize)/sigma)
                + 1 - tdist.cdf(tAlphaUpper + meanDiff*Math.sqrt(sampleSize)/sigma));
        
        return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
    }
    
    /*
     * Run a simulation of power for a one-tailed Student's t test.  The null and alternative distributions
     * are assumed to have equal variance.
     */
    protected OneSampleStudentsTPower simulateOneTailedPower(double mu0, double muA, 
    		double alpha, double sigma, int sampleSize, int simulationIterations)
    {       
        // calculate degrees of freedom
        int df = sampleSize-1;
        
        // create a central t distribution 
        StudentsT tdist = new StudentsT(df);
        // number of times the null hypothesis is correctly rejected during the simulation
        int rejections = 0;

        // run the simulation
        for(int i = 0; i < simulationIterations; i++)
        {
            // calculate a random value from the distribution under which the alternate hypothesis is true
            double x = sigma*tdist.random()/Math.sqrt(sampleSize) + muA;
            // convert to the corresponding t value in the distribution under which the alternative hypothesis is true
            double t0 = (x - mu0)*Math.sqrt(sampleSize)/sigma;
            // calculate the p value for t0 given that the null hypothesis is true
            double p = (t0 > 0 ? 1-tdist.cdf(t0) : tdist.cdf(t0));
            // accept or reject the null
            if (p <= alpha) rejections++;
        }

        double power = (double)rejections/(double)simulationIterations;
        
        return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
    }

    /*
     * Run a simulation of power for a two-tailed Student's t test.  The null and alternative distributions
     * are assumed to have equal variance.
     */
    protected OneSampleStudentsTPower simulateTwoTailedPower(double mu0, double muA, 
    		double alpha, double sigma, int sampleSize, int simulationIterations)
    {
        // calculate degrees of freedom
        int df = sampleSize-1;
        
        // create a t distribution from which we can draw random values 
        StudentsT tdist = new StudentsT(df);
        // number of times the null hypothesis is correctly rejected during the simulation
        int rejections = 0;

        // run the simulation
        for(int i = 0; i < simulationIterations; i++)
        {
            // calculate a random value from the distribution under which the null hypothesis is true
            double x = sigma*tdist.random()/Math.sqrt(sampleSize) + muA;
            // convert to the corresponding t value in the distribution under which the alternative hypothesis is true
            double t0 = (x - mu0)*Math.sqrt(sampleSize)/sigma;
            // calculate the p value for t0 given that the null hypothesis is true
            double p = (t0 > 0 ? 2*(1-tdist.cdf(t0)) : 2*(tdist.cdf(t0)));
            // accept or reject the null
            if (p <= alpha) rejections++;
        }

        double power = (double)rejections/(double)simulationIterations;
        
        return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
    }
    
    /**
     * Calculate sample size using a one tailed test
     * 
     * @param mu0 mean under the null hypothesis
     * @param muA mean under the alternative hypothesis
     * @param sigma standard deviation (assumed equal)
     * @param alpha type I error level
     * @param power desired power
     * @return sample size
     * @throws InvalidAlgorithmParameterException
     */
    private OneSampleStudentsTPower calculateOneTailedSampleSize(double mu0, double muA, double alpha,
            double sigma, double power) throws IllegalArgumentException
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(mu0, muA, sigma, alpha, power, true);
        
        try
        {
            int sampleSize = (int) Math.ceil(solver.solve(sampleSizeFunc, MIN_SAMPLE_SIZE, MAX_SAMPLE_SIZE));
            return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to calculate sample size: " + e.getMessage());
        }

    }

    /**
     * Calculate sample size using a two tailed test
     * 
     * @param mu0 mean under the null hypothesis
     * @param muA mean under the alternative hypothesis
     * @param sigma standard deviation (assumed equal)
     * @param alpha type I error level
     * @param power desired power
     * @return sample size
     * @throws InvalidAlgorithmParameterException
     */
    public OneSampleStudentsTPower calculateTwoTailedSampleSize(double mu0, double muA, double alpha,
            double sigma, double power) throws IllegalArgumentException
    {
        /* validate the inputs */
        if (sigma < 0)
            throw new IllegalArgumentException("Invalid standard error: " + sigma);
        if (power < 0 || power > 1)
            throw new IllegalArgumentException("Invalid power: " + power);
        if (alpha < 0 || alpha > 0.5)
            throw new IllegalArgumentException("Invalid alpha level: " + alpha);

        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();
        
        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(mu0, muA, sigma, alpha, power, false);
        
        try
        {
            int sampleSize = (int) Math.ceil(solver.solve(sampleSizeFunc, MIN_SAMPLE_SIZE, MAX_SAMPLE_SIZE));
            return new OneSampleStudentsTPower(alpha, power, sampleSize, Math.abs(mu0 - muA), sigma);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to calculate sample size: " + e.getMessage());
        }
    }
    
}
