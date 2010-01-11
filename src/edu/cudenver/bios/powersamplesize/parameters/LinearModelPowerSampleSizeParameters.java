package edu.cudenver.bios.powersamplesize.parameters;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;

/**
 * Container class for matrix inputs for general linear model power calculations.
 * <ul>
 * <li>Design matrix - essence matrix for the model design
 * <li>Beta matrix - estimated regression coefficients matrix
 * <li>Sigma matrix - estimated errors matrix
 * <li>Theta0 matrix - estimated null hypothesis matrix 
 * <li>Between subject contrast matrix - defines comparisons between subjects
 * <li>Within subject contrast matrix - defines comparisons within subjects
 * (may be left null for univariate special case)
 * </ul>
 * @author kreidles
 *
 */
public class LinearModelPowerSampleSizeParameters extends PowerSampleSizeParameters
{
    public enum TestStatistic {
        NONE,
        UNIREP,
        WILKS_LAMBDA,
        PILLAU_BARTLETT_TRACE,
        HOTELLING_LAWLEY_TRACE
    };
    
    public enum RandomPredictorAdjustment
    {
        CONDITIONAL_POWER,
        QUANTILE_POWER,
        INTEGRATION
    };
    
    TestStatistic testStatistic = TestStatistic.NONE;
    
    RealMatrix beta = null;
    
    RealMatrix betaOriginal = null;
    
    // used if only fixed predictors
    RealMatrix sigma = null;
    
    /* variance/ covariances required for baseline covariate */
    RealMatrix correlationCovariateOutcome = null;
    RealMatrix sigmaCovariate = null;
    RealMatrix sigmaOutcome = null;
    
    RealMatrix betweenSubjectContrast = null;
    
    RealMatrix theta = null;
    
    RealMatrix withinSubjectContrast = null;

    RealMatrix design = null;
    
    EssenceMatrix designEssence = null;
    
    RandomPredictorAdjustment randomPredictorAdjustment = 
        RandomPredictorAdjustment.CONDITIONAL_POWER;
    
    /**
     * Constructor.  Creates an empty set of linear model power parameters
     */
    public LinearModelPowerSampleSizeParameters() {}
    
    /**
     * Copy constructor. Creates a copy of the specified power parameters
     * 
     * @param params
     */
    public LinearModelPowerSampleSizeParameters(LinearModelPowerSampleSizeParameters params)
    {
        super(params);
        this.beta = params.getBeta();
        this.betaOriginal = new Array2DRowRealMatrix(beta.getData());
        this.design = params.getDesign();
        this.designEssence = params.getDesignEssence();
        this.sigma = params.getSigma();
        this.theta = params.getTheta();
        this.betweenSubjectContrast = params.getBetweenSubjectContrast();
        this.withinSubjectContrast = params.getWithinSubjectContrast();
        this.testStatistic = params.getTestStatistic();
        this.randomPredictorAdjustment = params.getRandomPredictorAdjustment();
    }
    
    public TestStatistic getTestStatistic()
    {
        return testStatistic;
    }

    public void setTestStatistic(TestStatistic testStatistic)
    {
        this.testStatistic = testStatistic;
    }

    public RealMatrix getBeta()
    {
        return beta;
    }

    public void setBeta(RealMatrix beta)
    {
        this.beta = beta;
        this.betaOriginal = new Array2DRowRealMatrix(beta.getData());
    }

    public RealMatrix getSigma()
    {
        return sigma;
    }

    public void setSigma(RealMatrix sigma)
    {
        this.sigma = sigma;
    }

    public RealMatrix getBetweenSubjectContrast()
    {
        return betweenSubjectContrast;
    }

    public void setBetweenSubjectContrast(RealMatrix betweenSubjectContrast)
    {
        this.betweenSubjectContrast = betweenSubjectContrast;
    }

    public RealMatrix getTheta()
    {
        return theta;
    }

    public void setTheta(RealMatrix theta)
    {
        this.theta = theta;
    }

    public RealMatrix getWithinSubjectContrast()
    {
        return withinSubjectContrast;
    }

    public void setWithinSubjectContrast(RealMatrix withinSubjectContrast)
    {
        this.withinSubjectContrast = withinSubjectContrast;
    }

    public RealMatrix getDesign()
    {
        if (design == null && designEssence != null)
        {
            // if only a design matrix is specified, retrieve the full design matrix
            this.design = designEssence.getFullDesignMatrix();
            // TODO: what to do if this is too big to cache?
        }
        return design;
    }
    
    public void setDesign(RealMatrix design)
    {
        this.design = design;
    }

    public EssenceMatrix getDesignEssence()
    {
        return designEssence;
    }

    public void setDesignEssence(EssenceMatrix designEssence)
    {
        this.designEssence = designEssence;
    }

    public RandomPredictorAdjustment getRandomPredictorAdjustment()
    {
        return randomPredictorAdjustment;
    }

    public void setRandomPredictorAdjustment(
            RandomPredictorAdjustment randomPredictorAdjustment)
    {
        this.randomPredictorAdjustment = randomPredictorAdjustment;
    }
    
    /**
     * Return the estimated sample size used in the power calculation
     * 
     * @return sample size
     */
    public int getSampleSize()
    {
    	int sampleSize = 0;
    	if (design != null)
    	{
    		sampleSize = design.getRowDimension();
    	}
    	else if (designEssence != null)
    	{
    		for (int i = 0; i < designEssence.getRowDimension(); i++)
    		{
        		RowMetaData rmd = designEssence.getRowMetaData(i);
    			sampleSize += rmd.getRepetitions();
    		}
    	}
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
        if (designEssence == null)
        	throw new IllegalArgumentException("Must specify an essence matrix to reset sample size");
        
        RealMatrix newDesign = designEssence.getFullDesignMatrix(sampleSize);
        design = newDesign;
    }
    
    /**
     * TODO: ask deb how to do this???
     */
    public void setMeanDifference(double delta) 
    throws IllegalArgumentException
    {
    	if (beta != null)
    	{
    		for(int row = 0; row < beta.getRowDimension(); row++)
    		{
    			for(int col = 0; col < beta.getColumnDimension(); col++)
    			{
    				if (betaOriginal.getEntry(row, col) != 0)  beta.setEntry(row, col, delta);
    			}
    		}
    	}
    }
    
    public double getMeanDifference()
    {
    	double delta = Double.NaN;
    	if (beta != null)
    	{
    		// TODO: wtf?
    	}
    	return delta;
    }
    
}
