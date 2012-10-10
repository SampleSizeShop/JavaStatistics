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
package edu.cudenver.bios.power.glmm;

import jsc.distributions.ChiSquared;
import jsc.distributions.FishersF;
import edu.cudenver.bios.distribution.NonCentralFDistribution;
import edu.cudenver.bios.power.PowerErrorEnum;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.utils.ConfidenceInterval;

public class GLMMPowerConfidenceInterval extends ConfidenceInterval
{
	public enum ConfidenceIntervalType
	{
		NONE,
		BETA_KNOWN_SIGMA_ESTIMATED,
		BETA_SIGMA_ESTIMATED
	}
	
	/**
	 * Compute confidence limits for power the GLMM(F) 
	 * (i.e. fixed predictors only).  This function is only called if the 
	 * ConfidenceInterval field in the parameters is not set to NONE
	 * 
	 * Based on the methods of Taylor and Muller:
	 * Taylor DJ and Muller KE (1995). Computing confidence bounds for power and sample
	 * size for the general linear univariate model.  The American Statistician. 49(1): 43-47;
	 * 
	 *  
	 * @param params
	 * @return list of confidence limits: lower limit, upper limit 
	 */
	public GLMMPowerConfidenceInterval(
			ConfidenceIntervalType ciType, double alphaLower, double alphaUpper, 
			int sampleSizeForEstimates, int designRankForEstimates, double alpha,
			GLMMTest test)
	throws PowerException
	{	
		// bail if no confidence limits are requested
		if (ciType == ConfidenceIntervalType.NONE)
			throw new PowerException("invalid confidence interval type",
			        PowerErrorEnum.POWER_CI_UNKNOWN_TYPE);
		// TODO: error checking
		
		double powerLower = alpha;
		double powerUpper = 1;

		if (test instanceof GLMMTestUnivariateRepeatedMeasures && test.isMultivariate())
		{
			// special case for the unirep CI's 
			// TODO: polymorphosize this somehow rather using a giant if/else
			// I don't think the theory exists for estimated beta and sigma for power CI's
			if (ciType == ConfidenceIntervalType.BETA_SIGMA_ESTIMATED && test.isMultivariate())
				throw new PowerException("cannot compute confidence limits for both beta and " + 
						"sigma estimated as the theory has not yet been derived",
						PowerErrorEnum.POWER_CI_MULTIVARIATE_BETA_SIGMA_ESTIMATED);

			/*
			 * Based on theoretical results for UNIREP from:
			 * Gribbin MJ (2007). Better Power Methods for the Univariate Approach to Repeated Measures.
			 * Ph.D. thesis, University of North Carolina at Chapel Hill.
			 * (see p65-68, equations 78-81)
			 */
			GLMMTestUnivariateRepeatedMeasures unirepTest = 
				(GLMMTestUnivariateRepeatedMeasures) test;
			
			// get the critical F under the null
			double FCrit = unirepTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
			// build the chi sqaure distribution for confidence limits
			double confLimitsDF = unirepTest.getConfidenceLimitsDegreesOfFreedom();
			ChiSquared chiSquaredDist = new ChiSquared(confLimitsDF);
			// get the lower bound for power
			if (alphaLower > 0)
			{
				double chiSquaredLower = chiSquaredDist.inverseCdf(alphaLower);
				double noncentralityLower = 
					(unirepTest.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE) * 
							chiSquaredLower / confLimitsDF);
				NonCentralFDistribution noncentralF = 
					new NonCentralFDistribution(unirepTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE), 
							unirepTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE), 
							noncentralityLower);
				powerLower = 1 - noncentralF.cdf(FCrit);
				//if (powerLower < alpha) powerLower = alpha;
				if (powerLower > 1) powerLower = 1;
			}

			if (alphaUpper > 0)
			{
				double chiSquaredUpper = chiSquaredDist.inverseCdf(1-alphaUpper);
				double noncentralityUpper = 
					(unirepTest.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE) * 
							chiSquaredUpper / confLimitsDF);
				NonCentralFDistribution noncentralF = 
					new NonCentralFDistribution(unirepTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE), 
							unirepTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE), 
							noncentralityUpper);
				powerUpper = 1 - noncentralF.cdf(FCrit);
				//if (powerUpper < alpha) powerUpper = alpha;
				if (powerUpper > 1) powerUpper = 1;
			}	
		}
		else
		{
			// get the degrees of freedom and noncentrality for the GLH we are using for power
			// (called the "target" sample)
			double targetHypothesisDF = test.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
			double targetErrorDF = test.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
			double omega = test.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE);
			double alternativeF = omega / targetHypothesisDF;

			// calculate the noncentrality for the estimation sample
			double estimationErrorDF = sampleSizeForEstimates - designRankForEstimates;
			double estimationNoncentrality = omega;
			FishersF centralFDist = new FishersF(targetHypothesisDF, targetErrorDF);
			double criticalF = centralFDist.inverseCdf(1 - alpha);

			// calculate the lower bound on the noncentrality 
			double noncentralityLower = 0;
			if (alphaLower > 0)
			{
				if (ciType == 	ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED)
				{
					ChiSquared chiSquareDist = new ChiSquared(estimationErrorDF);
					double chiLower = chiSquareDist.inverseCdf(alphaLower);
					noncentralityLower = (chiLower / estimationErrorDF) * estimationNoncentrality;
				}
				else if (ciType == ConfidenceIntervalType.BETA_SIGMA_ESTIMATED)
				{
					FishersF boundFDist = new FishersF(targetHypothesisDF, estimationErrorDF);
					double boundLower = boundFDist.inverseCdf(1 - alphaLower);
					if (alternativeF <= boundLower)
					{
						noncentralityLower = 0;
					}
					else
					{
						try
						{
							noncentralityLower = NonCentralFDistribution.noncentrality(alternativeF, 1- alphaLower,
									targetHypothesisDF, estimationErrorDF, estimationNoncentrality);
						}
						catch (Exception e)
						{
							// TODO better handling here?
							noncentralityLower = 0;
						}
					}
				}
			}

			// calculate the lower bound for power

			if (alphaLower > 0)
			{
				NonCentralFDistribution powerFDist = 
					new NonCentralFDistribution(targetHypothesisDF, targetErrorDF, noncentralityLower);
				powerLower = 1- powerFDist.cdf(criticalF);
			}
			if (powerLower < alpha)
			{
				// minimum power is the alpha level
				powerLower = alpha;
			}

			// now let's work on the upper limit
			// calculate the upper bound for the noncentrality
			double noncentralityUpper = Double.POSITIVE_INFINITY;
			if (alphaUpper > 0)
			{
				if (ciType == ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED)
				{
					ChiSquared chiSquareDist = new ChiSquared(estimationErrorDF);
					double chiUpper = chiSquareDist.inverseCdf(1 - alphaUpper);
					noncentralityUpper = (chiUpper / estimationErrorDF) * estimationNoncentrality;
				}
				else if (ciType == ConfidenceIntervalType.BETA_SIGMA_ESTIMATED)
				{
					FishersF boundFDist = new FishersF(targetHypothesisDF, estimationErrorDF);
					double boundUpper = boundFDist.inverseCdf(alphaUpper);
					if (alternativeF <= boundUpper)
					{
						noncentralityUpper = 0;
					}
					else
					{
						try
						{
							noncentralityUpper = NonCentralFDistribution.noncentrality(alternativeF, alphaUpper,
									targetHypothesisDF, estimationErrorDF, estimationNoncentrality);
						}
						catch (Exception e)
						{
							// TODO better handling here?
							noncentralityLower = 0;
						}
					}
				}
			}

			// calculate the upper power bound

			if (alphaUpper > 0)
			{
				NonCentralFDistribution powerFDist = 
					new NonCentralFDistribution(targetHypothesisDF, targetErrorDF, noncentralityUpper);
				powerUpper = 1- powerFDist.cdf(criticalF);
			}
			if (powerUpper < alpha)
			{
				// minimum power is the alpha level
				powerUpper = alpha;
			}
		}
		
		// store the confidence limit information
		this.alphaLower = alphaLower;
		this.alphaUpper = alphaUpper;
		this.upperLimit = powerUpper;
		this.lowerLimit = powerLower;

	}

}
