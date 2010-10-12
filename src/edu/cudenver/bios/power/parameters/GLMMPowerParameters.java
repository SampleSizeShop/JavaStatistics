/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2010 Regents of the University of Colorado.  
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package edu.cudenver.bios.power.parameters;

import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;

/**
 * Container class for inputs for general linear model power calculations.
 * A power/sample size calculation required matrix inputs to define the study
 * design, and lists defining the tests, type I error levels, scale factors, etc to
 * include in the power calculations.  
 * <p>
 * The following matrices should be specified
 * <ul>
 * <li>Design matrix - (required) essence matrix for the model design</li>
 * <li>Beta matrix - (required) estimated regression coefficients matrix.  When a baseline
 * covariate is specified, this matrix will have both fixed and random components</li>
 * <li>Theta matrix - (required) estimated null hypothesis matrix </li>
 * <li>Between subject contrast matrix - (required) defines comparisons between subjects. When a baseline
 * covariate is specified, this matrix will have both fixed and random components</li>
 * <li>Within subject contrast matrix - (required for multivariate) defines comparisons within subjects</li>
 * </ul>
 * <p></p>
 * For designs with fixed effects only, a single covariance matrix should be specified
 * <ul>
 * <li>Sigma (error)  matrix - estimated covariance matrix for the residual error</li>
 * </ul>
 * For designs with a baseline covariate, three covariance matrices should be specified:
 * <ul>
 * <li>Sigma Outcome (Y) - covariance matrix of responses</li>
 * <li>Sigma Gaussian Random (G) - covariance matrix of the baseline covariate (1x1)</li>
 * <li>SIgma Outcome Gaussian Random (YG) - covariance between responses and baseline covariate</li>
 * </ul>
 * The user may perform multiple power calculations on the above set of
 * matrices by specifying lists (with the corresponding "add" function) of the following values:
 * <ul>
 * <li>alpha - type I error values</li>
 * <li>group sample size - size of each group</li>
 * <li>power - desired power (for sample size calculations)</li>
 * <li>test - statistical test (ex. Wilk's Lambda, Hotelling-Lawley Trace, etc.)</li>
 * <li>betaScale - scale factors for the beta matrix.  Allows tweaking of mean difference estimates</li>
 * <li>sigmaScale - scale factors for the sigma error matrix.  Allows tweaking of residual covariance estimates</li>
 * <li>Power method - for designs with a baseline covariate, the user request unconditional or quantile power</li>
 * <li>Quantile - for quantile power calculations, the user may enter various quantiles (ex. 0.5 = median power)</li>
 * </ul>
 * @see PowerParameters
 * @see edu.cudenver.bios.matrix.DesignEssenceMatrix
 * @see edu.cudenver.bios.matrix.FixedRandomMatrix
 * @author Sarah Kreidler
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

	// power methods
	public enum PowerMethod
	{
		CONDITIONAL_POWER,
		UNCONDITIONAL_POWER,
		QUANTILE_POWER
	};

	// type of approximation to use for unirep
	public enum UnivariateCdf
	{
		MULLER_BARTON_APPROX,
		MULLER_EDWARDS_TAYLOR_APPROX,
		MULLER_EDWARDS_TAYLOR_EXACT,
		MULLER_EDWARDS_TAYLOR_EXACT_APPROX
	};

	// available moment approximation methods
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
	
	// residual error matrix and associated scale factors
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

	// set internally, full design matrix
	RealMatrix design = null;
	// caching of X'X inverse and rank since these are order(n^3)
	RealMatrix XtXInverse = null;
	int designRank = -1;
	// the essence design matrix
	DesignEssenceMatrix designEssence = null;

	// for design matrices with a baseline covariate, power may be estimated
	// with either conditional power (same as fixed effects), quantile power,
	// or unconditional power
	PeekableList<PowerMethod> powerMethodList = new PeekableList<PowerMethod>();
	// list of quantiles to use for quantile power
	PeekableList<Double> quantileList = new PeekableList<Double>();

	// default approximation method
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

	/**
	 * Add the statistical test to those included in the power/sample size 
	 * calculations
	 * @param test statistical test
	 */
	public void addTest(Test test)
	{
    	testList.add(test);
	}

	/**
	 * Lookup the default approximation for the current test
	 * @param test statistical test
	 */
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
	
	/**
	 * Internally activates the first statistical test and prepares
	 * for iteration over the test list
	 * @return first test
	 */
    public Test getFirstTest()
    {
    	Test test = testList.first();
    	setTestDefaults(test);
        return test;
    }
    
	/**
	 * Iterates to and activates the next  statistical test in the list.
	 * @return next test or null if at end of list
	 */
    public Test getNextTest()
    {
    	Test test = testList.next();
    	setTestDefaults(test);
        return test;
    }
    
    /**
     * Peek at the currently active statistical test
     * @return currently active test
     */
    public Test getCurrentTest()
    {
        return testList.current();
    }
	
    /**
     * Get the beta matrix
     * @return beta matrix
     */
	public FixedRandomMatrix getBeta()
	{
		return beta;
	}

	/**
	 * Set the estimated regression coefficients matrix
	 * @param beta beta matrix
	 */
	public void setBeta(FixedRandomMatrix beta)
	{
		this.beta = beta;
	}

	/**
	 * Get the residual error covariance matrix
	 * @return sigma error
	 */
	public RealMatrix getSigmaError()
	{
		return sigmaError;
	}

	/**
	 * Set the residual error covariance matrix
	 * @param sigmaError residual error covariance matrix
	 */
	public void setSigmaError(RealMatrix sigmaError)
	{
		this.sigmaError = sigmaError;
	}

	/**
	 * Get covariance matrix of responses and baseline covariate
	 * @return sigma Y-G
	 */
	public RealMatrix getSigmaOutcomeGaussianRandom()
	{
		return sigmaOutcomeGaussianRandom;
	}

	/**
	 * Set the covariance matrix of responses and baseline covariate
	 * @param sigmaOutcomeGaussianRandom covariance Y-G matrix
	 */
	public void setSigmaOutcomeGaussianRandom(RealMatrix sigmaOutcomeGaussianRandom)
	{
		this.sigmaOutcomeGaussianRandom = sigmaOutcomeGaussianRandom;
	}

	/**
	 * Get the covariance matrix of the baseline covariate (1x1)
	 * @return sigma G
	 */
	public RealMatrix getSigmaGaussianRandom()
	{
		return sigmaGaussianRandom;
	}

	/**
	 * Set the covariance matrix of the baseline covariate (1x1)
	 * @param sigmaGaussianRandom sigma G
	 */
	public void setSigmaGaussianRandom(RealMatrix sigmaGaussianRandom)
	{
		this.sigmaGaussianRandom = sigmaGaussianRandom;
	}

	/**
	 * Get the covariance of the responses
	 * @return sigma Y
	 */
	public RealMatrix getSigmaOutcome()
	{
		return sigmaOutcome;
	}

	/**
	 * Set the covariance of the responses
	 * @param sigmaOutcome sigma Y
	 */
	public void setSigmaOutcome(RealMatrix sigmaOutcome)
	{
		this.sigmaOutcome = sigmaOutcome;
	}

	/**
	 * Get the between subject contrast matrix (C)
	 * @return C matrix
	 */
	public FixedRandomMatrix getBetweenSubjectContrast()
	{
		return betweenSubjectContrast;
	}

	/**
	 * Set the between subject contrast matrix
	 * @param betweenSubjectContrast C matrix
	 */
	public void setBetweenSubjectContrast(FixedRandomMatrix betweenSubjectContrast)
	{
		this.betweenSubjectContrast = betweenSubjectContrast;
	}

	/**
	 * Get the theta null matrix
	 * @return theta null
	 */
	public RealMatrix getTheta()
	{
		return theta;
	}

	/**
	 * Set the theta null matrix
	 * @param theta null hypothesis matrix
	 */
	public void setTheta(RealMatrix theta)
	{
		this.theta = theta;
	}

	/**
	 * Get the within subject contrast (U)
	 * @return U matrix
	 */
	public RealMatrix getWithinSubjectContrast()
	{
		return withinSubjectContrast;
	}

	/**
	 * Set the within subject contrast (U)
	 * @param withinSubjectContrast U matrix
	 */
	public void setWithinSubjectContrast(RealMatrix withinSubjectContrast)
	{
		this.withinSubjectContrast = withinSubjectContrast;
	}

	/**
	 * Get the full design matrix.  This value is cached.
	 * @return design matrix
	 */
	public RealMatrix getDesign()
	{
		return getDesign(false);
	}
	
	/**
	 * Get rank of design matrix
	 * 
	 * @return rank of design matrix
	 */
	public int getDesignRank()
	{
		if (designRank < 0)
		{
			designRank = new SingularValueDecompositionImpl(getDesign()).getRank();
		}
		return designRank;
	}
	
	/**
	 * Get cached (X'X)-1.  Matrix inversion is order n^3 so we
	 * cache this value
	 * 
	 * @return (X'X) inverse
	 */
	public RealMatrix getXtXInverse()
	{
		if (XtXInverse == null)
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
			XtXInverse = null;
			designRank = -1;
		}
		else if (force);
		{
			designEssence.updateRandomColumns(design);
			XtXInverse = null;
		}
		return design;
	}
	
	/**
	 * Set the full design matrix - in most cases, this is set
	 * internally and an essence matrix should be used.
	 * @param design
	 */
	public void setDesign(RealMatrix design)
	{
		this.design = design;
	}

	/**
	 * Get the "essence" design matrix which contains 1 entry for
	 * each row in the fixed portion of the full design matrix 
	 * @return essence design matrix
	 */
	public DesignEssenceMatrix getDesignEssence()
	{
		return designEssence;
	}

	/**
	 * Set the "essence" design matrix which contains 1 entry for
	 * each row in the fixed portion of the full design matrix 
	 * @param designEssence design essence matrix.
	 */
	public void setDesignEssence(DesignEssenceMatrix designEssence)
	{
		this.designEssence = designEssence;
	}

	/**
	 * Get the current moment approximation method
	 * @return moment approximation method
	 */
	public MomentApproximationMethod getMomentMethod()
	{
		return momentMethod;
	}

	/**
	 * Set the current moment approximation method
	 * @param momentMethod moment approximation method
	 */
	public void setMomentMethod(MomentApproximationMethod momentMethod)
	{
		this.momentMethod = momentMethod;
	}

	/**
	 * Get the method for estimating the univariate cdf
	 * @return univariate cdf method
	 */
	public UnivariateCdf getUnivariateCdf()
	{
		return univariateCdf;
	}

	/**
	 * Set the method for estimating the univariate cdf
	 * @param univariateCdf unvariate cdf method
	 */
	public void setUnivariateCdf(UnivariateCdf univariateCdf)
	{
		this.univariateCdf = univariateCdf;
	}

    /**
     * Add a scale factor for the sigma error matrix.  Allows the user
     * to tweak the estimates of residual variance in the model
     * @param scale scale factor
     */
    public void addSigmaScale(double scale)
    {
    	sigmaScaleList.add(new Double(scale));    	
    }
    
    /**
     * Add a scale factor for the beta matrix.  This scale factor
     * only effects the fixed portion of the beta matrix
     * @param scale scale factor
     */
    public void addBetaScale(double scale)
    {
    	betaScaleList.add(new Double(scale));    	
    }    
	
	/**
	 * Internally scales the beta matrix with the first beta scale and prepares
	 * for iteration over the beta scale list
	 * @return first beta scale factor
	 */
    public Double getFirstBetaScale()
    {
        Double betaScale = betaScaleList.first();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue(), true);
        return betaScale;
    }
    
	/**
	 * Internally scales the beta matrix with the next beta scale in the beta 
	 * scale list
	 * @return next beta scale factor
	 */
    public Double getNextBetaScale()
    {
        Double betaScale = betaScaleList.next();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue(), true);
        return betaScale;
    }
    
	/**
	 * Peek at the current beta scale factor
	 * @return current beta scale factor
	 */
    public Double getCurrentBetaScale()
    {
        return betaScaleList.current();
    }
	
    /**
     * Get the scaled beta matrix.  Value depends on previous calls to
     * getFirstBetaScale and getNextBetaScale
     * @return scaled beta matrix
     */
	public RealMatrix getScaledBeta()
	{
	    return (betaScaled != null ? betaScaled : beta.getCombinedMatrix());
	}
	
	/**
	 * Internally scales the sigma error matrix with the first sigma scale and prepares
	 * for iteration over the sigma scale list
	 * @return first sigma scale factor
	 */
    public Double getFirstSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.first();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;        
    }
	/**
	 * Internally scales the sigma matrix with the next sigma scale in the sigma 
	 * scale list
	 * @return next sigma scale factor
	 */
    public Double getNextSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.next();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;    
    }
	/**
	 * Peek at the current sigma scale factor
	 * @return current sigma scale factor
	 */
    public Double getCurrentSigmaScale()
    {
        return sigmaScaleList.current();
    }
    
    /**
     * Get the scaled sigma matrix.  Value depends on previous calls to
     * getFirstSigmaScale and getNextSigmaScale
     * @return scaled sigma matrix
     */
    public RealMatrix getScaledSigmaError()
    {
        return (sigmaErrorScaled != null ? sigmaErrorScaled : sigmaError);
    }
    
	/**
	 * Internally selects the first per group sample size and sets the 
	 * design to null, forcing the full design matrix to be regenerated when
	 * next requested
	 * @return first per group sample size
	 */
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
    
	/**
	 * Internally selects the next per group sample size and sets the 
	 * design to null, forcing the full design matrix to be regenerated when
	 * next requested
	 * @return next per group sample size
	 */
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
    
    /**
     * Peek at the current sample size
     */
    @Override
    public Integer getCurrentSampleSize()
    {
        return sampleSizeList.current();
    }
    
    /**
     * Set the per-group sample size
     * @param sampleSize
     */
    public void setGroupSampleSize(int sampleSize)
    {
        designEssence.setGroupSampleSize(sampleSize);
        design = null;
    }

    /**
     * Get the list of beta scale factors
     * @return beta scale list
     */
    public PeekableList<Double> getBetaScaleList()
    {
        return betaScaleList;
    }

    /**
     * Get the list of sigma scale factors
     * @return sigma scale list
     */
    public PeekableList<Double> getSigmaScaleList()
    {
        return sigmaScaleList;
    }

    /**
     * Force a rescale of the beta matrix
     * @param betaScale scale factor
     */
    public void scaleBeta(double betaScale)
    {
        if (betaScale != Double.NaN)
            betaScaled = beta.scalarMultiply(betaScale, true);
    }
    
    /**
     * Reset the scaled beta matrix
     * @param scaledBeta new scaled beta matrix
     */
    public void setScaledBeta(RealMatrix scaledBeta)
    {
        betaScaled = scaledBeta;
    }
    
    /**
     * Force a rescale of the sigma error matrix
     * @param sigmaScale scale factor
     */
    public void scaleSigmaError(double sigmaScale)
    {
        if (sigmaScale != Double.NaN)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale);
    }
    
    /**
     * Reset the scaled sigma error matrix
     * @param scaledSigma new scaled sigma error matrix
     */
    public void setScaledSigmaError(RealMatrix scaledSigma)
    {
        sigmaErrorScaled = scaledSigma;
    }
    
    /**
     * Add a power method for use in power calculations
     * @param method conditional, quantile, or unconditional power
     */
    public void addPowerMethod(PowerMethod method)
    {
    	powerMethodList.add(method);
    }
    
	/**
	 * Internally selects the first power method and prepares the list
	 * for iteration
	 * @return first power method
	 */
    public PowerMethod getFirstPowerMethod()
    {
    	PowerMethod method = powerMethodList.first();
    	// if no methods are specified, just use conditional power
    	if (method == null)
    		return DEFAULT_POWER_METHOD;
    	else
    		return method;
    }
    
    /**
     * Activates the next power method 
     * @return next power method
     */
    public PowerMethod getNextPowerMethod()
    {
    	PowerMethod method = powerMethodList.next();
        return method;
    }
    
    /**
     * Peeks at the currently active power method
     * @return currently active power method
     */
    public PowerMethod getCurrentPowerMethod()
    {
    	PowerMethod method = powerMethodList.current();
    	// if no methods are specified, just use conditional power
    	if (method == null)
    		return DEFAULT_POWER_METHOD;
    	else
    		return method;
    		
    }
        
    /**
     * Add a quantile (value between 0 and 1) to the list for use with the quantile
     * power method.
     * @param quantile value between 0 and 1 indicating the power quantile
     */
    public void addQuantile(double quantile)
    {
    	quantileList.add(new Double(quantile));
    }
    
	/**
	 * Internally selects the first quantile and prepares the list
	 * for iteration.  Only used for quantile power.
	 * @return first quantile
	 */
    public Double getFirstQuantile()
    {
    	Double quantile = quantileList.first();
    	// if no quantiles specified, return default
    	if (quantile == null)
    		return new Double(DEFAULT_QUANTILE); 
    	else
    		return quantile;
    }
    
    /**
     * Activates the next quantile in the quantile list
     * @return next quantile
     */
    public Double getNextQuantile()
    {
    	Double quantile = quantileList.next();
        return quantile;
    }
    
    /**
     * Peek at currently active quantile value
     * @return current quantile
     */
    public Double getCurrentQuantile()
    {
        return quantileList.current();
    }

    /**
     * Indicates if the estimate of the non-centrality cdf
     * uses the exact Davies' algorithm.  If false, a Satterthwaite
     * approximation is used.
     * @return true if exact Davies's algorithm is used.
     */
	public boolean isNonCentralityCDFExact() 
	{
		return nonCentralityCDFExact;
	}

    /**
     * Set whether the estimate of the non-centrality cdf
     * should use the exact Davies' algorithm.  If false, a Satterthwaite
     * approximation is used.
     * 
     * @param nonCentralityCDFExact true for Davies' algorithm, false
     * for Satterthwaite approximation.
     */
	public void setNonCentralityCDFExact(boolean nonCentralityCDFExact) 
	{
		this.nonCentralityCDFExact = nonCentralityCDFExact;
	}
}


