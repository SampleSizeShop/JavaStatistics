package edu.cudenver.bios.power.parameters;

import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;

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
public class GLMMPowerParameters extends PowerParameters
{
	private static final PowerMethod DEFAULT_POWER_METHOD = 
		PowerMethod.CONDITIONAL_POWER;
	private static final double DEFAULT_QUANTILE = 0.5;
	
	// the type of statistical test to use
	public enum Test 
	{
		UNIREP,
		UNIREP_BOX,
		UNIREP_GEISSER_GREENHOUSE,
		UNIREP_HUYNH_FELDT,
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

	// type of statistical test being performed
	PeekableList<Test> testList = new PeekableList<Test>();
	
	// beta matrix of regression coefficients
	FixedRandomMatrix beta = null;
	RealMatrix betaScaled = null; // beta matrix with scale factor applied
	PeekableList<Double> betaScaleList = new PeekableList<Double>();
	
	// used if only fixed predictors
	RealMatrix sigmaError = null;
	RealMatrix sigmaErrorScaled = null;
	PeekableList<Double> sigmaScaleList = new PeekableList<Double>();
	
	/* variance/ covariances required for random gaussian covariates */
	RealMatrix sigmaOutcomeGaussianRandom = null;
	RealMatrix sigmaGaussianRandom = null;
	RealMatrix sigmaOutcome = null;

	// C matrix - contrasts for between subject effects
	FixedRandomMatrix betweenSubjectContrast = null;

	// theta matrix - matrix of null hypothesis values
	RealMatrix theta = null;

	// U matrix - contrasts for within subject effects
	RealMatrix withinSubjectContrast = null;

	RealMatrix design = null;
	RealMatrix XtXInverse = null;
	int designRank = -1;
	DesignEssenceMatrix designEssence = null;

	// for design matrices with a baseline covariate, power may be estimated
	// with either conditional power (same as fixed effects), quantile power,
	// or unconditional power
	PeekableList<PowerMethod> powerMethodList = new PeekableList<PowerMethod>();
	// list of quantiles to use for quantile power
	PeekableList<Double> quantileList = new PeekableList<Double>();

	UnivariateCdf univariateCdf = UnivariateCdf.MULLER_EDWARDS_TAYLOR_APPROX;

	MomentApproximationMethod momentMethod =
		MomentApproximationMethod.NONE;

	// if true, use the exact calculation of the CDF of the non-centrality parameter
	// (applies to quantile and conditional power only)
	boolean nonCentralityCDFExact = false;
	
	/**
	 * Constructor.  Creates an empty set of linear model power parameters
	 */
	public GLMMPowerParameters() 
	{
	    super();
	}

	public void addTest(Test test)
	{
    	testList.add(test);
	}

	private void setTestDefaults(Test test)
	{
		if (test != null)
		{
			switch (test)
			{
			case HOTELLING_LAWLEY_TRACE:
				momentMethod = MomentApproximationMethod.MCKEON_TWO_MOMENT_OMEGA_MULT;
				break;
			case PILLAI_BARTLETT_TRACE:
				momentMethod = MomentApproximationMethod.MULLER_TWO_MOMENT;
				break;
			case WILKS_LAMBDA:
				momentMethod = MomentApproximationMethod.RAO_TWO_MOMENT_OMEGA_MULT;
				break;
			default:
				momentMethod = MomentApproximationMethod.NONE;
			}
		}
	}
	
    public Test getFirstTest()
    {
    	Test test = testList.first();
    	setTestDefaults(test);
        return test;
    }
    
    public Test getNextTest()
    {
    	Test test = testList.next();
    	setTestDefaults(test);
        return test;
    }
    
    public Test getCurrentTest()
    {
        return testList.current();
    }
	
	public FixedRandomMatrix getBeta()
	{
		return beta;
	}

	public void setBeta(FixedRandomMatrix beta)
	{
		this.beta = beta;
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

	public FixedRandomMatrix getBetweenSubjectContrast()
	{
		return betweenSubjectContrast;
	}

	public void setBetweenSubjectContrast(FixedRandomMatrix betweenSubjectContrast)
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
			// if only an essence matrix is specified, retrieve the full design matrix
			this.design = designEssence.getFullDesignMatrix();
			// TODO: what to do if this is too big to cache?
		}
		return design;
	}
	
	/**
	 * Added caching of rank of design matrix
	 * 
	 * @param force - if true, force recomputation of (X'X)-1
	 * @return
	 */
	public int getDesignRank(boolean force)
	{
		if (designRank < 0 || force)
		{
			designRank = new SingularValueDecompositionImpl(getDesign()).getRank();
		}
		return designRank;
	}
	
	/**
	 * Added caching of (X'X)-1 since this is order n^3
	 * 
	 * @param force - if true, force recomputation of (X'X)-1
	 * @return
	 */
	public RealMatrix getXtXInverse(boolean force)
	{
		if (XtXInverse == null || force)
		{
			RealMatrix X = getDesign();
	        XtXInverse = new LUDecompositionImpl(X.transpose().multiply(X)).getSolver().getInverse();
		}
		return XtXInverse;
	}
	
	/**
	 * Regenerates the design matrix and fills any random columns with a new
	 * realization of random values based on a normal distribution.
	 * 
	 * @return full design matrix
	 */
	public RealMatrix getDesign(boolean force)
	{
		if (design == null && designEssence != null)
		{
			// if only an essence matrix is specified, retrieve the full design matrix
			design = designEssence.getFullDesignMatrix();
		}
		else
		{
			designEssence.updateRandomColumns(design);
		}
		return design;
	}
	
	public void setDesign(RealMatrix design)
	{
		this.design = design;
	}

	public DesignEssenceMatrix getDesignEssence()
	{
		return designEssence;
	}

	public void setDesignEssence(DesignEssenceMatrix designEssence)
	{
		this.designEssence = designEssence;
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

    /**
     * Add a scale factor for the sigma matrix
     */
    public void addSigmaScale(double scale)
    {
    	sigmaScaleList.add(new Double(scale));    	
    }
    
    /**
     * Add a scale factor for the beta matrix
     */
    public void addBetaScale(double scale)
    {
    	betaScaleList.add(new Double(scale));    	
    }    
	
    public Double getFirstBetaScale()
    {
        Double betaScale = betaScaleList.first();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue(), true);
        return betaScale;
    }
    
    public Double getNextBetaScale()
    {
        Double betaScale = betaScaleList.next();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue(), true);
        return betaScale;
    }
    
    public Double getCurrentBetaScale()
    {
        return betaScaleList.current();
    }
	
	public RealMatrix getScaledBeta()
	{
	    return (betaScaled != null ? betaScaled : beta.getCombinedMatrix());
	}
	
	
    public Double getFirstSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.first();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;        
    }
    
    public Double getNextSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.next();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;    
    }
    
    public Double getCurrentSigmaScale()
    {
        return sigmaScaleList.current();
    }
    
    public RealMatrix getScaledSigmaError()
    {
        return (sigmaErrorScaled != null ? sigmaErrorScaled : sigmaError);
    }
    
    @Override
    public Integer getFirstSampleSize()
    {
        Integer sampleSize = sampleSizeList.first();
        if (sampleSize != null)
        {
            designEssence.setGroupSampleSize(sampleSize.intValue());
            design = null;
        }
        return sampleSize;
    }
    
    @Override
    public Integer getNextSampleSize()
    {
        Integer sampleSize = sampleSizeList.next();
        if (sampleSize != null)
        {
            designEssence.setGroupSampleSize(sampleSize.intValue());
            design = null;
        }
        return sampleSize;
    }
    
    @Override
    public Integer getCurrentSampleSize()
    {
        return sampleSizeList.current();
    }
    
    public void setGroupSampleSize(int sampleSize)
    {
        designEssence.setGroupSampleSize(sampleSize);
        design = null;
    }

    public PeekableList<Double> getBetaScaleList()
    {
        return betaScaleList;
    }

    public RealMatrix getSigmaErrorScaled()
    {
        return sigmaErrorScaled;
    }

    public PeekableList<Double> getSigmaScaleList()
    {
        return sigmaScaleList;
    }

    public void scaleBeta(double betaScale)
    {
        if (betaScale != Double.NaN)
            betaScaled = beta.scalarMultiply(betaScale, true);
    }
    
    public void setScaledBeta(RealMatrix scaledBeta)
    {
        betaScaled = scaledBeta;
    }
    
    public void scaleSigmaError(double sigmaScale)
    {
        if (sigmaScale != Double.NaN)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale);
    }
    
    public void setScaledSigmaError(RealMatrix scaledSigma)
    {
        sigmaErrorScaled = scaledSigma;
    }
    
    public void addPowerMethod(PowerMethod method)
    {
    	powerMethodList.add(method);
    }
    
    public PowerMethod getFirstPowerMethod()
    {
    	PowerMethod method = powerMethodList.first();
    	// if no methods are specified, just use conditional power
    	if (method == null)
    		return DEFAULT_POWER_METHOD;
    	else
    		return method;
    }
    
    public PowerMethod getNextPowerMethod()
    {
    	PowerMethod method = powerMethodList.next();
        return method;
    }
    
    public PowerMethod getCurrentPowerMethod()
    {
    	PowerMethod method = powerMethodList.current();
    	// if no methods are specified, just use conditional power
    	if (method == null)
    		return DEFAULT_POWER_METHOD;
    	else
    		return method;
    		
    }
        
    public void addQuantile(double quantile)
    {
    	quantileList.add(new Double(quantile));
    }
    
    public Double getFirstQuantile()
    {
    	Double quantile = quantileList.first();
    	// if no quantiles specified, return default
    	if (quantile == null)
    		return new Double(DEFAULT_QUANTILE); 
    	else
    		return quantile;
    }
    
    public Double getNextQuantile()
    {
    	Double quantile = quantileList.next();
        return quantile;
    }
    
    public Double getCurrentQuantile()
    {
        return quantileList.current();
    }

	public boolean isNonCentralityCDFExact() 
	{
		return nonCentralityCDFExact;
	}

	public void setNonCentralityCDFExact(boolean nonCentralityCDFExact) 
	{
		this.nonCentralityCDFExact = nonCentralityCDFExact;
	}
}


