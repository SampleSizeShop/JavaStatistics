package edu.cudenver.bios.matrix;

/**
 * Container class for column specific meta data in an essence matrix.
 * Allows the user to specify  a mean and variance of the predictor 
 * distribution for a given column
 * 
 * @author Sarah Kreidler
 *
 */
public class RandomColumnMetaData
{
	double mean = Double.NaN;
    double variance = Double.NaN;
    
    /**
     * Constructor.  Create a column meta data object with the specified predictor
     * type.  Mean and variance should be specified for random predictors.
     * 
     * @param predictorType indicates if the predictor is fixed or random
     * @param mean mean of the predictor distribution (assumed normal)
     * @param variance variance of the predictor distribution
     */
    public RandomColumnMetaData (double mean, double variance) 
    {
        this.mean = mean;
        this.variance = variance;
    }

    /**
     * Get the mean for this column meta data type.  This value is
     * valid only for random predictors.
     * 
     * @return mean of the predictor distribution
     */
    public double getMean()
    {
        return mean;
    }

    /**
     * Set the mean for the predictor distribution.  Only valid 
     * for random predictors.
     * 
     * @param mean mean of the predictor distribution
     */
    public void setMean(double mean)
    {
        this.mean = mean;
    }

    /**
     * Get the variance of the predictor distribution.  Only valid 
     * for random predictors.
     * 
     * @return variance of predictor distribution
     */
    public double getVariance()
    {
        return variance;
    }

    /**
     * Set the variance of the predictor distribution.  Only valid for
     * random predictors.
     * 
     * @param variance variance of the predictor distribution
     */
    public void setVariance(double variance)
    {
        this.variance = variance;
    }
    
}
