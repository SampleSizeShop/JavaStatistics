package edu.cudenver.bios.powersamplesize;

import java.security.InvalidAlgorithmParameterException;

import jsc.distributions.StudentsT;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;

import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.SimplePowerSampleSizeParameters;

/**
 * Class to estimate sample size for a one sample Student's T test
 * 
 * @author Sarah Kreidler
 *
 */
public class SampleSizeOneSampleStudentsT implements SampleSize
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
    
    /**
     * Constructor
     */
    public SampleSizeOneSampleStudentsT() {}
    
    /**
     * Estimate the sample size for the Student's T test for the given parameters
     * 
     * @see SimpleSampleSizeParameters
     * 
     * @param params SimpleSampleSizeParameters object
     */
    public int getSampleSize(PowerSampleSizeParameters params)
    throws IllegalArgumentException
    {
        SimplePowerSampleSizeParameters ssp = (SimplePowerSampleSizeParameters) params;
        
        // validate inputs
        validateParameters(ssp);

        // calculate the power for either one or two tails
        if (ssp.isOneTailed())
        {
            return calculateOneTailedSampleSize(ssp.getMu0(), ssp.getMuA(), ssp.getAlpha(), ssp.getSigma(), ssp.getPower());
        }
        else
        {
            return calculateTwoTailedSampleSize(ssp.getMu0(), ssp.getMuA(), ssp.getAlpha(), ssp.getSigma(), ssp.getPower());
        }
    }
    
    /*
     * Validate the sample size parameters
     */
    private void validateParameters(SimplePowerSampleSizeParameters ssp)
    throws IllegalArgumentException
    {
        if (!ssp.isInitializedMu0()|| !ssp.isInitializedMuA())
            throw new IllegalArgumentException("Must specify both null and alternative means");
        if (ssp.getSigma() < 0)
            throw new IllegalArgumentException("Invalid standard error: " + ssp.getSigma());
        if (ssp.getPower() <= 0 || ssp.getPower() >= 1)
            throw new IllegalArgumentException("Invalid pwoer: " + ssp.getPower());
        if (ssp.getAlpha() < 0 || ssp.getAlpha() > 1)
            throw new IllegalArgumentException("Invalid alpha level: " + ssp.getAlpha());
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
    private int calculateOneTailedSampleSize(double mu0, double muA, double alpha,
            double sigma, double power) throws IllegalArgumentException
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(mu0, muA, sigma, alpha, power, true);
        
        try
        {
            return (int) Math.ceil(solver.solve(sampleSizeFunc, MIN_SAMPLE_SIZE, MAX_SAMPLE_SIZE));
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
    public int calculateTwoTailedSampleSize(double mu0, double muA, double alpha,
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
            return (int) Math.ceil(solver.solve(sampleSizeFunc, MIN_SAMPLE_SIZE, MAX_SAMPLE_SIZE));
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to calculate sample size: " + e.getMessage());
        }
    }
}
