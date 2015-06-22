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

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateEpsilonApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;

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
	// power methods
	public enum PowerMethod
	{
		CONDITIONAL_POWER,
		UNCONDITIONAL_POWER,
		QUANTILE_POWER
	};

	// parameters related to confidence intervals
	ConfidenceIntervalType confidenceIntervalType = ConfidenceIntervalType.NONE;
	double alphaLowerConfidenceLimit = 0.025;
	double alphaUpperConfidenceLimit = 0.025;
	int sampleSizeForEstimates = 0;
	int designMatrixRankForEstimates = 0;

	/** flags for selecting approximation methods, cdf methods **/
	// approximation method setting.  These can be set for each statistical test
	HashMap<Test,UnivariateCdfApproximation> univariateCdfMap =
		new HashMap<Test,UnivariateCdfApproximation>();
	HashMap<Test,UnivariateEpsilonApproximation> univariateEpsilonMap =
		new HashMap<Test,UnivariateEpsilonApproximation>();
	HashMap<Test,FApproximation> FMap =
		new HashMap<Test,FApproximation>();
	// if true, use the exact calculation of the CDF of the non-centrality parameter
	// (applies to quantile and conditional power only)
	boolean nonCentralityCDFExact = false;

	/* lists which define a set of powers calculations for the current run
	 * note that sample size, power, and alpha lists are defined in the PowerParameters super class
	 */
	// type of statistical test being performed
	ArrayList<Test> testList = new ArrayList<Test>();
	// scale factors for the beta matrix
	ArrayList<Double> betaScaleList = new ArrayList<Double>();
	// scale factors for the sigma(error) matrix
	ArrayList<Double> sigmaScaleList = new ArrayList<Double>();
	// for design matrices with a baseline covariate, power may be estimated
	// with either conditional power (same as fixed effects), quantile power,
	// or unconditional power
	ArrayList<PowerMethod> powerMethodList = new ArrayList<PowerMethod>();
	// list of quantiles to use for quantile power
	ArrayList<Double> quantileList = new ArrayList<Double>();

	/* matrix inputs to the power calculation.  */
	// the design essence matrix
	/* For details please see Muller & Fetterman (2002) "Regression and ANOVA" */
	RealMatrix designEssence = null;
	// caching of X'X inverse and rank since these are order(n^3) operations
	RealMatrix XtXInverse = null;
	int designRank = -1;
	// C matrix - contrasts for between subject effects
	FixedRandomMatrix betweenSubjectContrast = null;
	// U matrix - contrasts for within subject effects
	RealMatrix withinSubjectContrast = null;
	// theta matrix - matrix of null hypothesis values
	RealMatrix theta = null;
	// beta matrix of regression coefficients
	FixedRandomMatrix beta = null;
	// residual error matrix
	RealMatrix sigmaError = null;
	// variance/ covariances required for random gaussian covariates
	// i.e. GLMM(F,g) designs
	RealMatrix sigmaOutcomeGaussianRandom = null;
	RealMatrix sigmaGaussianRandom = null;
	RealMatrix sigmaOutcome = null;

	/**
	 * Constructor.  Creates an empty set of linear model power parameters
	 */
	public GLMMPowerParameters()
	{
	    super();
	    setTestDefaults();
	}

	/**** Functions for list inputs which define multiple power calculations ****/

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
     * Add a scale factor for the beta matrix.  This scale factor
     * only effects the fixed portion of the beta matrix
     * @param scale scale factor
     */
    public void addBetaScale(double scale)
    {
    	betaScaleList.add(new Double(scale));
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
     * Add a power method for use in power calculations
     * @param method conditional, quantile, or unconditional power
     */
    public void addPowerMethod(PowerMethod method)
    {
    	powerMethodList.add(method);
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
	 * Get the list of statistical tests
	 * @return list of statistical tests requested
	 */
	public ArrayList<Test> getTestList()
	{
		return testList;
	}

    /**
     * Get the list of beta scale factors
     * @return beta scale list
     */
    public ArrayList<Double> getBetaScaleList()
    {
        return betaScaleList;
    }

    /**
     * Get the list of sigma scale factors
     * @return sigma scale list
     */
    public ArrayList<Double> getSigmaScaleList()
    {
        return sigmaScaleList;
    }

    /**
     * Get the list of power methods
     * @return power method list
     */
    public ArrayList<PowerMethod> getPowerMethodList()
    {
    	return powerMethodList;
    }

    /**
     * Get the list of quantiles associated with a set of quantile power calculations
     * @return list of quantiles
     */
    public ArrayList<Double> getQuantileList()
    {
    	return quantileList;
    }

	/**
	 * Clear the list of statistical tests
	 */
	public void clearTestList()
	{
		testList.clear();
	}

    /**
     * Clear the list of beta scale factors
     */
    public void clearBetaScaleList()
    {
        betaScaleList.clear();
    }

    /**
     * Clear the list of sigma scale factors
     */
    public void clearSigmaScaleList()
    {
        sigmaScaleList.clear();
    }

    /**
     * Clear the list of power methods
     */
    public void clearPowerMethodList()
    {
    	powerMethodList.clear();
    }

    /**
     * Clear the list of quantiles associated with a set of quantile power calculations
     */
    public void clearQuantileList()
    {
    	quantileList.clear();
    }

    /**
     * Clear the list of per group sample size values
     */
    public void clearSampleSizeList()
    {
        sampleSizeList.clear();
    }

	/**** Functions for setting approximation and cdf methods ****/

    /**
     * Sets the approximation method for the F statistic which is created
     * by the underlying statistical test.
     * @param test the statistical test for which the F approximation is being set
     */
    public void setFApproximationMethod(Test test, FApproximation method)
    {
    	FMap.put(test, method);
    }

    /**
     * Get the F approximation method for the specified test
     * @param test
     */
    public FApproximation getFApproximationMethod(Test test)
    {
    	return FMap.get(test);
    }

    /**
     * Sets the CDF approximation method for the univariate tests
     * @param test statistical test type
     * @param method CDF approximation method
     */
    public void setUnivariateCdfMethod(Test test, UnivariateCdfApproximation method)
    {
    	univariateCdfMap.put(test, method);
    }

    /**
     * Get the CDF approximation method for the specified test
     * @param test statistical test
     * @return CDF approximation method
     */
    public UnivariateCdfApproximation getUnivariateCdfMethod(Test test)
    {
    	return univariateCdfMap.get(test);
    }

    /**
     * Sets the approximation method for mean epsilon the univariate tests
     * @param test statistical test type
     * @param method epsilon approximation method
     */
    public void setUnivariateEpsilonMethod(Test test, UnivariateEpsilonApproximation method)
    {
    	univariateEpsilonMap.put(test, method);
    }

    /**
     * Get the approximation method for mean epsilon for the specified univariate test
     * @param test statistical test
     * @return Epsilon approximation method
     */
    public UnivariateEpsilonApproximation getUnivariateEpsilonMethod(Test test)
    {
    	return univariateEpsilonMap.get(test);
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

	/**** Functions for setting the matrix inputs ****/

	/**
	 * Get the "essence" design matrix which contains 1 entry for
	 * each row in the fixed portion of the full design matrix
	 * @return essence design matrix
	 */
	public RealMatrix getDesignEssence()
	{
		return designEssence;
	}

	/**
	 * Set the "essence" design matrix which contains 1 entry for
	 * each row in the fixed portion of the full design matrix
	 * @param designEssence design essence matrix.
	 */
	public void setDesignEssence(RealMatrix designEssence)
	{
		this.designEssence = designEssence;
		designRank = new SingularValueDecomposition(designEssence).getRank();
		if (this.sigmaGaussianRandom != null) designRank++;
		XtXInverse =
        	new LUDecomposition(this.designEssence.transpose().multiply(this.designEssence)).getSolver().getInverse();
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
			designRank = new SingularValueDecomposition(designEssence).getRank();
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
	        XtXInverse =
	        	new LUDecomposition(designEssence.transpose().multiply(designEssence)).getSolver().getInverse();
		}
		return XtXInverse;
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
		if (this.sigmaGaussianRandom == null) designRank++;
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
	 * Fill in the F approximation method defaults
	 */
	private void setTestDefaults()
	{
		for(GLMMTestFactory.Test test : GLMMTestFactory.Test.values())
		{
			switch (test)
			{
			case HOTELLING_LAWLEY_TRACE:
				FMap.put(test, FApproximation.MCKEON_TWO_MOMENT_OMEGA_MULT);
				break;
			case PILLAI_BARTLETT_TRACE:
				FMap.put(test, FApproximation.MULLER_TWO_MOMENT);
				break;
			case WILKS_LAMBDA:
				FMap.put(test, FApproximation.RAO_TWO_MOMENT_OMEGA_MULT);
				break;
			default:
				FMap.put(test, FApproximation.NONE);
			}
		}
	}

    /**** Confidence interval settings ****/

	/**
	 * Get the lower tail probability on the power confidence interval
	 */
	public double getAlphaLowerConfidenceLimit()
	{
		return alphaLowerConfidenceLimit;
	}

	/**
	 * Set the lower tail probability for the power confidence interval
	 * @param alphaLowerConfidenceLimit lower tail probability for the confidence interval
	 */
	public void setAlphaLowerConfidenceLimit(double alphaLowerConfidenceLimit)
	{
		this.alphaLowerConfidenceLimit = alphaLowerConfidenceLimit;
	}

	/**
	 * Get the upper tail probability on the power confidence interval
	 */
	public double getAlphaUpperConfidenceLimit()
	{
		return alphaUpperConfidenceLimit;
	}

	/**
	 * Set the upper tail probability for the power confidence interval
	 * @param alphaUpperConfidenceLimit upper tail probability for the confidence interval
	 */
	public void setAlphaUpperConfidenceLimit(double alphaUpperConfidenceLimit)
	{
		this.alphaUpperConfidenceLimit = alphaUpperConfidenceLimit;
	}

	/**
	 * Get the sample size for the data set from which the beta and
	 * sigma estimates were derived
	 * @return sample size for data set used to estimate beta and sigma
	 */
	public int getSampleSizeForEstimates()
	{
		return sampleSizeForEstimates;
	}

	/**
	 * Set the sample size of the data set from which the beta and sigma estimates
	 * were derived
	 * @param sampleSizeForEstimates sample size for data set used to estimate beta and sigma
	 */
	public void setSampleSizeForEstimates(int sampleSizeForEstimates)
	{
		this.sampleSizeForEstimates = sampleSizeForEstimates;
	}

	/**
	 * Get the rank of the design matrix for the data set used to estimate beta and sigma
	 * @return rank of the design matrix for the estimation data set
	 */
	public int getDesignMatrixRankForEstimates()
	{
		return designMatrixRankForEstimates;
	}

	/**
	 * Set the rank of the design matrix for the data set used to estimate beta and sigma
	 * @param designMatrixRankForEstimates rank of the design matrix for the estimation data set
	 */
	public void setDesignMatrixRankForEstimates(int designMatrixRankForEstimates)
	{
		this.designMatrixRankForEstimates = designMatrixRankForEstimates;
	}

	/**
	 * Get the type of confidence interval.  The confidence intervals can either be constructed
	 * under the assumption that both beta and sigma are estimates,
	 * or that beta is fixed and only sigma is estimated.
	 * @return type of confidence interval
	 */
	public ConfidenceIntervalType getConfidenceIntervalType()
	{
		return confidenceIntervalType;
	}

	/**
	 * Get the type of confidence interval.  The confidence intervals can either be constructed
	 * under the assumption that both beta and sigma are estimates,
	 * or that beta is fixed and only sigma is estimated.
	 * @param ConfidenceIntervalType type of confidence interval
	 */
	public void setConfidenceIntervalType(ConfidenceIntervalType ConfidenceIntervalType)
	{
		this.confidenceIntervalType = ConfidenceIntervalType;
	}
}
