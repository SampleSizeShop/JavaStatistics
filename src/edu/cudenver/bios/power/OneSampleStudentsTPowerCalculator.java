package edu.cudenver.bios.power;

import java.security.InvalidAlgorithmParameterException;
import java.util.ArrayList;
import java.util.List;

import jsc.distributions.StudentsT;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;

import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters.MeanPair;

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
        boolean twoTailed = false;
        
        public SampleSizeFunction(double mu0, double muA, double sigma,
                double alpha, double power, boolean twoTailed)
        {
            this.mu0 = mu0;
            this.muA = muA;
            this.sigma = sigma;
            this.alpha = alpha;
            this.power = power;
            this.twoTailed = twoTailed;
        }
        
        public double value(double n)
        {
            try
            {
                // create a t distribution with the specified degrees of freedom
                double df = n -1;
                StudentsT tdist = new StudentsT(df);

                if (twoTailed)
                {
                    double tAlpha = tdist.inverseCdf(1-alpha/2);
                    double tBeta = tdist.inverseCdf(power);
                    double root = ((sigma*sigma*(tBeta + tAlpha)*(tBeta + tAlpha))/
                            ((mu0 - muA)*(mu0-muA))) - n;
                    return root;
                }
                else
                {
                    double tAlpha = tdist.inverseCdf(1-alpha);
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
        for(Double alpha = studentsTParams.getFirstAlpha(); alpha != null;
                alpha = studentsTParams.getNextAlpha())
        {
            for(Double sigma = studentsTParams.getFirstSigma(); sigma != null;
                    sigma = studentsTParams.getNextSigma())
            {
                for(MeanPair means = studentsTParams.getFirstMeans(); means != null;
                        means = studentsTParams.getNextMeans())
                {
                    for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
                            sampleSize = params.getNextSampleSize())
                    {
                        try
                        {
                            results.add(calculatePower(alpha.doubleValue(), means.mu0, means.muA, 
                                    sigma.doubleValue(), sampleSize.intValue(), studentsTParams.isTwoTailed()));
                        }
                        catch (MathException me)
                        {
                            // TODO: what to do?
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

	    // calculate the sample size for either one or two tails
	    for(Double alpha = studentsTParams.getFirstAlpha(); alpha != null;
	    alpha = studentsTParams.getNextAlpha())
	    {
	        for(Double sigma = studentsTParams.getFirstSigma(); sigma != null;
	        sigma = studentsTParams.getNextSigma())
	        {
	            for(MeanPair means = studentsTParams.getFirstMeans(); means != null;
	            means = studentsTParams.getNextMeans())
	            {
	                for(Double power = params.getFirstPower(); power != null; 
	                power = params.getNextPower())
	                {
	                    results.add(calculateSampleSize(alpha.doubleValue(), means.mu0, means.muA, 
	                            sigma.doubleValue(), power.doubleValue(), studentsTParams.isTwoTailed()));
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
        for(Double alpha = studentsTParams.getFirstAlpha(); alpha != null;
        alpha = studentsTParams.getNextAlpha())
        {
            for(Double sigma = studentsTParams.getFirstSigma(); sigma != null;
            sigma = studentsTParams.getNextSigma())
            {
                for(MeanPair means = studentsTParams.getFirstMeans(); means != null;
                means = studentsTParams.getNextMeans())
                {
                    for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
                    sampleSize = params.getNextSampleSize())
                    {
        				try
        				{
        					results.add(simulatePower(alpha.doubleValue(), means.mu0, means.muA,
        							sigma.doubleValue(), sampleSize.intValue(), studentsTParams.isTwoTailed(), iterations));
        				}
        				catch (MathException me)
        				{
        					// TODO: how to handle this?
        				}
        			}
        		}
        	}
        }

        return results;
	}

	protected OneSampleStudentsTPower calculatePower(double alpha, double mu0, double muA, 
			double sigma, int sampleSize, boolean isTwoTail) throws MathException
	{
        // calculate the degrees of freedom 
        int df = sampleSize - 1;
        // get the critical t for specified alpha level in H0 distribution (central T) 
        StudentsT tdist = new StudentsT(df);
    	// calculate the difference of means 
        double meanDiff = Math.abs(mu0 - muA);
        // the resulting power value
        double power = Double.NaN;
        
        if (isTwoTail)
        {
            double tAlphaLower = tdist.inverseCdf(alpha/2);
            double tAlphaUpper = tdist.inverseCdf(1-alpha/2);
            power = (tdist.cdf(tAlphaLower + meanDiff*Math.sqrt(sampleSize)/sigma)
                    + 1 - tdist.cdf(tAlphaUpper + meanDiff*Math.sqrt(sampleSize)/sigma));
        }
        else
        {        	
            double tAlpha = tdist.inverseCdf(alpha);
            power = (tdist.cdf(tAlpha + meanDiff*Math.sqrt(sampleSize)/sigma));
        }
        	
        return new OneSampleStudentsTPower(alpha, power, power, sampleSize, meanDiff, sigma);
	
	}
    
	 /**
     * Calculate sample size for the one sample t test.  This function uses a bisection search algorithm 
     * to determine sample size.  The actual power is also calculated.
     * 
     * @param mu0 mean under the null hypothesis
     * @param muA mean under the alternative hypothesis
     * @param sigma standard deviation (assumed equal)
     * @param alpha type I error level
     * @param power desired power
     * @return sample size
     * @throws InvalidAlgorithmParameterException
     */
	protected OneSampleStudentsTPower calculateSampleSize(double alpha, double mu0, double muA, 
    		double sigma, double power, boolean isTwoTailed)
	{
        /* validate the inputs */
        if (sigma < 0)
            throw new IllegalArgumentException("Invalid standard error: " + sigma);
        if (power < 0 || power > 1)
            throw new IllegalArgumentException("Invalid power: " + power);
        if (alpha <= 0 || alpha >= 1)
            throw new IllegalArgumentException("Invalid alpha level: " + alpha);
        
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(mu0, muA, sigma, alpha, power, isTwoTailed);
        
        try
        {
            int sampleSize = (int) Math.ceil(solver.solve(sampleSizeFunc, MIN_SAMPLE_SIZE, MAX_SAMPLE_SIZE));
            OneSampleStudentsTPower actualPower = calculatePower(alpha, mu0, muA, sigma, sampleSize, isTwoTailed);
            actualPower.setNominalPower(power);
            return actualPower;
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to calculate sample size: " + e.getMessage());
        }
	}
	
	protected OneSampleStudentsTPower simulatePower(double alpha, double mu0, double muA, 
    		double sigma, int sampleSize, boolean isTwoTailed, int simulationIterations) throws MathException
	{
        // calculate degrees of freedom
        int df = sampleSize-1;
        // create a central t distribution 
        StudentsT tdist = new StudentsT(df);

        // number of times the null hypothesis is correctly rejected during the simulation
        int rejections = 0;
        
        if (isTwoTailed)
        {
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
        }
        else 
        {
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
        }
        
        double power = (double)rejections/(double)simulationIterations;
        
        return new OneSampleStudentsTPower(alpha, power, power, sampleSize, Math.abs(mu0 - muA), sigma);
	}
    
}
