package edu.cudenver.bios.matrix;

/**
 * Container class for column specific meta data in an essence matrix.
 * Allows the user to specify if the column represents a fixed or random 
 * predictor.  If the predictor is random, the user may specify a mean and 
 * variance of the predictor distribution.
 * 
 * @author Sarah Kreidler
 *
 */
public class ColumnMetaData
{
    // types of predictor variables (columns of design matrix)
    public enum PredictorType 
    {
      FIXED,
      RANDOM
    };
    
    PredictorType predictorType = PredictorType.FIXED;
    double mean = Double.NaN;
    double variance = Double.NaN;

    /**
     * Constructor.  Defaults to fixed predictor type.
     */
    public ColumnMetaData () {}
    
    /**
     * Constructor.  Create a column meta data object with the specified predictor
     * type.  Mean and variance should be specified for random predictors.
     * 
     * @param predictorType indicates if the predictor is fixed or random
     * @param mean mean of the predictor distribution (assumed normal)
     * @param variance variance of the predictor distribution
     */
    public ColumnMetaData (PredictorType predictorType,
            double mean, double variance) 
    {
        this.predictorType = predictorType;
        this.mean = mean;
        this.variance = variance;
    }
    
    /**
     * Get the predictor type (fixed or random)
     * 
     * @return predictor type
     */
    public PredictorType getPredictorType()
    {
        return predictorType;
    }

    /**
     * Set the predictor type for this column meta data
     * 
     * @param predictorType
     */
    public void setPredictorType(PredictorType predictorType)
    {
        this.predictorType = predictorType;
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
