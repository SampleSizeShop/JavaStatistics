package edu.cudenver.bios.powersamplesize.parameters;

/**
 * Container class for power input parameters.  Most statistical tests
 * will need to extend this class.
 * 
 * @author Sarah Kreidler
 *
 */
public abstract class PowerSampleSizeParameters
{    
    /**
     * Desired power
     */
    double power = -1;
    
    /**
     * Alpha level (type I error level)
     */
    double alpha;


    
    /**
     * Create an empty power parameter object
     */
    public PowerSampleSizeParameters() {}

    /**
     * Copy constructor
     */
    public PowerSampleSizeParameters(PowerSampleSizeParameters params) 
    {
        this.alpha = params.getAlpha();
        this.power = params.getPower();
    }
    
    /**
     * Get alpha level
     * @return alpha level
     */
    public double getAlpha()
    {
        return alpha;
    }
    
    /**
     * Set alpha level
     * 
     * @param alpha
     */
    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
    }
    


    /**
     * Get the desired power 
     * 
     * @return power
     */
    public double getPower()
    {
        return power;
    }

    /**
     * Set the desired power
     * 
     * @param power
     */
    public void setPower(double power)
    {
        this.power = power;
    }
    
    /**
     * Set the mean difference.  
     * 
     * @param delta
     */
    abstract public void setMeanDifference(double delta) throws IllegalArgumentException;
    
    /**
     * Get the mean difference
     * 
     * @return mean difference
     */
    abstract public double getMeanDifference();
    
    /**
     * Get the current sample size
     * 
     * @return
     */
    abstract public int getSampleSize();
    
    /**
     * Set the sample size
     * 
     * @param n
     */
    abstract public void setSampleSize(int n) throws IllegalArgumentException;
}

