package edu.cudenver.bios.powersamplesize;

import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;

/**
 * Interface representing calculations and simulations of statistical power
 * 
 * @author Sarah Kreidler
 */
public interface Power
{
    /**
     * Calculate power for a statistical test
     * 
     * @param params collection of inputs to the power calculation
     * @return calculated power value
     * @throws IllegalArgumentException
     */
    public double getCalculatedPower(PowerSampleSizeParameters params) 
    throws IllegalArgumentException;
    
    /**
     * Runs a simulation of power for a statistical test with the
     * specified number of iterations.  Simulations run the following steps:
     * <ul>
     * <li>Generate a random value from the appropriate distribution for the statistical test
     * <li>Calculate the p-value for the given test
     * <li>Count the number of times the p-value is small enough to reject the null hypothesis
     * <li>The frequency of rejections is the simulated power
     * </ul>
     * 
     * @param params collection of inputs to the power simulation
     * @param iterations number of times to run a simulated
     * @return calculated power value
     * @throws IllegalArgumentException
     */
    public double getSimulatedPower(PowerSampleSizeParameters params, int iterations)
    throws IllegalArgumentException;
    
    
}
