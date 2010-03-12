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
        PILLAI_BARTLETT_TRACE,
        HOTELLING_LAWLEY_TRACE
    };
    
    public enum PowerMethod
    {
        CONDITIONAL_POWER,
        UNCONDITIONAL_POWER,
        QUANTILE_POWER
    };
    
    public enum UnivariateCorrection 
    {
        NONE,
        BOX,
        GEISSER_GREENHOUSE,
        HUYNH_FELDT
    };
    
    public enum UnivariateCdf
    {
        MULLER_BARTON_APPROX,
        MULLER_EDWARDS_TAYLOR_APPROX,
        MULLER_EDWARDS_TAYLOR_EXACT,
        MULLER_EDWARDS_TAYLOR_EXACT_APPROX
    };
    
    public enum MomentApproximationMethod
    {
        NONE,
        PILLAI_ONE_MOMENT,
        PILLAI_ONE_MOMENT_OMEGA_MULT,
        MCKEON_TWO_MOMENT,
        MCKEON_TWO_MOMENT_OMEGA_MULT,
        MULLER_TWO_MOMENT,
        MULLER_TWO_MOMENT_OMEGA_MULT,
        RAO_TWO_MOMENT,
        RAO_TWO_MOMENT_OMEGA_MULT
    };
    
    TestStatistic testStatistic = TestStatistic.NONE;
    
    RealMatrix beta = null;
    
    RealMatrix betaOriginal = null;
    
    // used if only fixed predictors
    RealMatrix sigmaError = null;
    
    /* variance/ covariances required for random gaussian covariates */
    RealMatrix sigmaOutcomeGaussianRandom = null;
    RealMatrix sigmaGaussianRandom = null;
    RealMatrix sigmaOutcome = null;
    
    RealMatrix betweenSubjectContrast = null;
    
    RealMatrix theta = null;
    
    RealMatrix withinSubjectContrast = null;

    RealMatrix design = null;
    
    EssenceMatrix designEssence = null;
    
    PowerMethod powerMethod = PowerMethod.CONDITIONAL_POWER;
    double quantile = 0.50;
    
    UnivariateCorrection univariateCorrection = UnivariateCorrection.NONE;
    
    UnivariateCdf univariateCdf = UnivariateCdf.MULLER_BARTON_APPROX;
    
    MomentApproximationMethod momentMethod =
        MomentApproximationMethod.NONE;
    
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
        this.setBeta(params.getBeta());
        this.design = params.getDesign();
        this.designEssence = params.getDesignEssence();
        this.sigmaError = params.getSigmaError();
        this.sigmaGaussianRandom = params.getSigmaGaussianRandom();
        this.sigmaOutcome = params.getSigmaOutcome();
        this.sigmaOutcomeGaussianRandom = params.getSigmaOutcomeGaussianRandom();
        this.theta = params.getTheta();
        this.betweenSubjectContrast = params.getBetweenSubjectContrast();
        this.withinSubjectContrast = params.getWithinSubjectContrast();
        this.testStatistic = params.getTestStatistic();
        this.powerMethod = params.getPowerMethod();
        this.quantile = params.getQuantile();
        this.univariateCorrection = params.getUnivariateCorrection();
        this.univariateCdf = params.getUnivariateCdf();
        this.momentMethod = params.getMomentMethod();
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
        if (beta != null)
            this.betaOriginal = new Array2DRowRealMatrix(beta.getData());
        else
            this.betaOriginal = null;
    }

    public RealMatrix getSigmaError()
    {
        return sigmaError;
    }

    public void setSigmaError(RealMatrix sigmaError)
    {
        this.sigmaError = sigmaError;
    }
    
    public RealMatrix getSigmaOutcomeGaussianRandom()
    {
        return sigmaOutcomeGaussianRandom;
    }

    public void setSigmaOutcomeGaussianRandom(RealMatrix sigmaOutcomeGaussianRandom)
    {
        this.sigmaOutcomeGaussianRandom = sigmaOutcomeGaussianRandom;
    }

    public RealMatrix getSigmaGaussianRandom()
    {
        return sigmaGaussianRandom;
    }

    public void setSigmaGaussianRandom(RealMatrix sigmaGaussianRandom)
    {
        this.sigmaGaussianRandom = sigmaGaussianRandom;
    }

    public RealMatrix getSigmaOutcome()
    {
        return sigmaOutcome;
    }

    public void setSigmaOutcome(RealMatrix sigmaOutcome)
    {
        this.sigmaOutcome = sigmaOutcome;
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

    public PowerMethod getPowerMethod()
    {
        return powerMethod;
    }

    public void setPowerMethod(
            PowerMethod powerMethod)
    {
        this.powerMethod = powerMethod;
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

    public double getQuantile()
    {
        return quantile;
    }

    public void setQuantile(double quantile)
    {
        this.quantile = quantile;
    }

    public UnivariateCorrection getUnivariateCorrection()
    {
        return univariateCorrection;
    }

    public void setUnivariateCorrection(UnivariateCorrection univariateCorrection)
    {
        this.univariateCorrection = univariateCorrection;
    }

    public MomentApproximationMethod getMomentMethod()
    {
        return momentMethod;
    }

    public void setMomentMethod(MomentApproximationMethod momentMethod)
    {
        this.momentMethod = momentMethod;
    }

    public UnivariateCdf getUnivariateCdf()
    {
        return univariateCdf;
    }

    public void setUnivariateCdf(UnivariateCdf univariateCdf)
    {
        this.univariateCdf = univariateCdf;
    }
    
}
