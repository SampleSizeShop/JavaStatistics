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
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.utils.ConfidenceInterval;

public class GLMMPowerConfidenceInterval extends ConfidenceInterval
{
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
	public GLMMPowerConfidenceInterval(GLMMPowerParameters params)
	throws IllegalArgumentException
	{
		// create a test
        GLMMTest test = GLMMTestFactory.createGLMMTest(params);
		
		// bail if no confidence limits are requested
		if (params.getConfidenceIntervalType() == 
			GLMMPowerParameters.ConfidenceIntervalType.NONE)
			throw new IllegalArgumentException("invalid confidence interval type");

		// full error checking performed in validateInputs routine

		// get the alpha tail probabilities
		double alphaLower = params.getAlphaLowerConfidenceLimit();
		double alphaUpper = params.getAlphaUpperConfidenceLimit();

		// get the degrees of freedom and noncentrality for the GLH we are using for power
		// (called the "target" sample)
		double targetHypothesisDF = test.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
		double targetErrorDF = test.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
		double omega = test.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE);
		double alternativeF = omega / targetHypothesisDF;

		// calculate the noncentrality for the estimation sample
		double estimationErrorDF = params.getSampleSizeForEstimates() - params.getDesignMatrixRankForEstimates();
		double estimationNoncentrality = omega;
		FishersF centralFDist = new FishersF(targetHypothesisDF, targetErrorDF);
		double criticalF = centralFDist.inverseCdf(1 - params.getCurrentAlpha());

		// calculate the lower bound on the noncentrality 
		double noncentralityLower = 0;
		if (alphaLower > 0)
		{
			if (params.getConfidenceIntervalType() == 
				GLMMPowerParameters.ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED)
			{
				ChiSquared chiSquareDist = new ChiSquared(estimationErrorDF);
				double chiLower = chiSquareDist.inverseCdf(alphaLower);
				noncentralityLower = (chiLower / estimationErrorDF) * estimationNoncentrality;
			}
			else if (params.getConfidenceIntervalType() == 
				GLMMPowerParameters.ConfidenceIntervalType.BETA_SIGMA_ESTIMATED)
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
		double powerLower = params.getCurrentAlpha();
		if (alphaLower > 0)
		{
			NonCentralFDistribution powerFDist = 
				new NonCentralFDistribution(targetHypothesisDF, targetErrorDF, noncentralityLower);
			powerLower = 1- powerFDist.cdf(criticalF);
		}
		if (powerLower < params.getCurrentAlpha())
		{
			// minimum power is the alpha level
			powerLower = params.getCurrentAlpha();
		}

		// now let's work on the upper limit
		// calculate the upper bound for the noncentrality
		double noncentralityUpper = Double.POSITIVE_INFINITY;
		if (alphaUpper > 0)
		{
			if (params.getConfidenceIntervalType() == 
				GLMMPowerParameters.ConfidenceIntervalType.BETA_KNOWN_SIGMA_ESTIMATED)
			{
				ChiSquared chiSquareDist = new ChiSquared(estimationErrorDF);
				double chiUpper = chiSquareDist.inverseCdf(1 - alphaUpper);
				noncentralityUpper = (chiUpper / estimationErrorDF) * estimationNoncentrality;
			}
			else if (params.getConfidenceIntervalType() == 
				GLMMPowerParameters.ConfidenceIntervalType.BETA_SIGMA_ESTIMATED)
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
		double powerUpper = 1;
		if (alphaUpper > 0)
		{
			NonCentralFDistribution powerFDist = 
				new NonCentralFDistribution(targetHypothesisDF, targetErrorDF, noncentralityUpper);
			powerUpper = 1- powerFDist.cdf(criticalF);
		}
		if (powerUpper < params.getCurrentAlpha())
		{
			// minimum power is the alpha level
			powerUpper = params.getCurrentAlpha();
		}

		this.alphaLower = alphaLower;
		this.alphaUpper = alphaUpper;
		this.upperLimit = powerUpper;
		this.lowerLimit = powerLower;
	}

}
