package edu.cudenver.bios.powersamplesize.parameters;

/**
 * Power/sample size parameter object for simple tests such as the 
 * Z, Student's T.
 * <br>
 * Contains null mean, alternative mean, standard deviation,
 * sample size (if calculating power) or power (if calculating sample size),
 * and alpha level
 * 
 * @author Sarah Kreidler
 *
 */
public class SimplePowerSampleSizeParameters extends PowerSampleSizeParameters
{
    // mean under the null hypothesis
    double mu0 = 0;
    boolean initializedMu0 = false;
    // mean under the alternative hypothesis
    double muA = 0;
    boolean initializedMuA = false;
    // estimated population std dev
    double sigma = -1;
    /**
     * Indicates if a one or two tailed power calculation should be performed
     */
    boolean oneTailed = false;
    
    /**
     * Estimated sample size 
     */
    int sampleSize = -1;
    
    /**
     * 
     * @return
     */
    public double getMu0()
    {
        return mu0;
    }
    
    public void setMu0(double mu0)
    {
        this.mu0 = mu0;
        this.initializedMu0 = true;
    }
    
    public boolean isInitializedMu0()
    {
        return initializedMu0;
    }

    public double getMuA()
    {
        return muA;
    }
    
    public void setMuA(double muA)
    {
        this.muA = muA;
        this.initializedMuA = true;
    }
    
    public boolean isInitializedMuA()
    {
        return initializedMuA;
    }
    
    public double getSigma()
    {
        return sigma;
    }
    
    public void setSigma(double sigma)
    {
        this.sigma = sigma;
    }

    /**
     * True if one-tailed, false if two-tailed power calculations will be performed
     * (default is two-tailed)
     * 
     * @return sample size
     */
    public boolean isOneTailed()
    {
        return oneTailed;
    }

    /**
     * Set to true for one-tailed, false for two-tailed power calculations
     * 
     * @return sample size
     */
    public void setOneTailed(boolean oneTailed)
    {
        this.oneTailed = oneTailed;
    }
    
    /**
     * Return the estimated sample size used in the power calculation
     * 
     * @return sample size
     */
    public int getSampleSize()
    {
        return sampleSize;
    }

    /**
     * Specify the estimated sample size for use in the power calculation
     * 
     * @return sample size
     */
    public void setSampleSize(int sampleSize)
    throws IllegalArgumentException
    {
        this.sampleSize = sampleSize;
    }
    
    public void setMeanDifference(double delta) 
    throws IllegalArgumentException
    {
    	if (!initializedMu0) mu0 = 0;
    	muA = mu0 + delta;
    }
    
    public double getMeanDifference()
    {
    	return muA - mu0;
    }
}
