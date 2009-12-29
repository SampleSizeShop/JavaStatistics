package edu.cudenver.bios.powersamplesize;

import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;

/**
 * Interface representing calculations and simulations of statistical power
 *  
 * @author Sarah Kreidler
 */
public interface SampleSize
{
    /**
     * Estimate sample size for a statistical test.  
     * 
     * @param params collection of inputs to the sample size calculation
     * @return estimated sample size
     * @throws IllegalArgumentException
     */
    public int getSampleSize(PowerSampleSizeParameters params) 
    throws IllegalArgumentException;
}
