/**
 * Calculate power for Student's T Distribution
 * 
 * Version: 1.0.0
 * 
 * Date: October 16, 2009
 * 
 * Copyright: ???
 * 
 * @author kreidles
 * 
 */
package edu.cudenver.bios.powersamplesize;

import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.SimplePowerSampleSizeParameters;
import jsc.distributions.StudentsT;

/**
 * Power implementation for the Student's T test.  Provides both one and two tailed tests
 * 
 * @author Sarah Kreidler
 *
 */
public class PowerOneSampleStudentsT implements Power
{
    /**
     * Constructor
     */
    public PowerOneSampleStudentsT() {}

    /**
     * Calculates the power of the Student's T test for the given parameters
     * 
     * @see SimplePowerSampleSizeParameters
     * 
     * @param params SimplePowerParameter object
     */
    public double getCalculatedPower(PowerSampleSizeParameters params)
    throws IllegalArgumentException
    {
        SimplePowerSampleSizeParameters spp = (SimplePowerSampleSizeParameters) params;
        
        // validate inputs
        validateParameters(spp);

        // calculate the power for either one or two tails
        if (spp.isOneTailed())
        {
            return calculateOneTailedPower(spp.getMu0(), spp.getMuA(), spp.getAlpha(), spp.getSigma(), spp.getSampleSize());
        }
        else
        {
            return calculateTwoTailedPower(spp.getMu0(), spp.getMuA(), spp.getAlpha(), spp.getSigma(), spp.getSampleSize());
        }
    }

    /**
     * Runs a power simulation of the Student's T test for the given parameters
     * 
     * @see SimplePowerParameters, Power
     * 
     * @param params SimplePowerParameter object
     * @param iterations number of simulation iterations
     * 
     * @return simulated power value (frequency of correct rejections of the null)
     */
    public double getSimulatedPower(PowerSampleSizeParameters params, int iterations)
    throws IllegalArgumentException
    {
        SimplePowerSampleSizeParameters spp = (SimplePowerSampleSizeParameters) params;
        
        // validate inputs
        validateParameters(spp);

        // calculate the power for either one or two tails
        if (spp.isOneTailed())
        {
            return simulateOneTailedPower(spp.getMu0(), spp.getMuA(), spp.getAlpha(), 
                    spp.getSigma(), spp.getSampleSize(), iterations);
        }
        else
        {
            return simulateTwoTailedPower(spp.getMu0(), spp.getMuA(), spp.getAlpha(), 
                    spp.getSigma(), spp.getSampleSize(), iterations);
        }
    } 
    
    /*
     * Calculate a one-tailed power
     */
    protected double calculateOneTailedPower(double mu0, double muA, double alpha,
            double sigma, int sampleSize) throws IllegalArgumentException
    {
        /* calculate the degrees of freedom */
        int df = sampleSize - 1;
        /* get the critical t for specified alpha level in H0 distribution (central T) */
        StudentsT tdist = new StudentsT(df);
        double tAlpha = tdist.inverseCdf(alpha);
        
        return(tdist.cdf(tAlpha + Math.abs(mu0 - muA)*Math.sqrt(sampleSize)/sigma));
    }

    /*
     * Calculate a two-tailed power
     */
    protected double calculateTwoTailedPower(double mu0, double muA, double alpha,
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
        return(tdist.cdf(tAlphaLower + meanDiff*Math.sqrt(sampleSize)/sigma)
                + 1 - tdist.cdf(tAlphaUpper + meanDiff*Math.sqrt(sampleSize)/sigma));
    }
    
    /*
     * Run a simulation of power for a one-tailed Student's t test.  The null and alternative distributions
     * are assumed to have equal variance.
     */
    protected double simulateOneTailedPower(double mu0, double muA, double alpha, double sigma, 
            int sampleSize, int simulationIterations)
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
            if (p < alpha) rejections++;
        }

        return (double)rejections/(double)simulationIterations;
    }

    /*
     * Run a simulation of power for a two-tailed Student's t test.  The null and alternative distributions
     * are assumed to have equal variance.
     */
    protected double simulateTwoTailedPower(double mu0, double muA, double alpha, 
            double sigma, int sampleSize, int simulationIterations)
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
            if (p < alpha) rejections++;
        }

        return (double)rejections/(double)simulationIterations;
    }
    /*
     * Make sure the parameters are valid for the power calculation/simulation
     */
    private void validateParameters(SimplePowerSampleSizeParameters spp) throws IllegalArgumentException
    {
        if (!spp.isInitializedMu0()|| !spp.isInitializedMuA())
            throw new IllegalArgumentException("Must specify both null and alternative means");
        if (spp.getSigma() < 0)
            throw new IllegalArgumentException("Invalid standard error: " + spp.getSigma());
        if (spp.getSampleSize() <= 1)
            throw new IllegalArgumentException("Invalid sample size: " + spp.getSampleSize());
        if (spp.getAlpha() < 0 || spp.getAlpha() > 1)
            throw new IllegalArgumentException("Invalid alpha level: " + spp.getAlpha());
    }
}
