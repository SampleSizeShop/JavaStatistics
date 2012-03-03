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
package edu.cudenver.bios.power;

import java.io.Serializable;

import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.utils.ConfidenceInterval;

/**
 * Pojo containing a description of the general linear model power result
 * 
 * @author Sarah Kreidler
 */
public class GLMMPower extends Power implements Serializable
{
    private static final long serialVersionUID = -6289570391909037726L;
    // scale factor for beta matrix
	double betaScale;
	// scale factor for the sigma error matrix
	double sigmaScale;
	// statistical test performed
	GLMMTestFactory.Test test;
	
	/* optional parameters */
	// power method - specified if a baseline covariate is used
	GLMMPowerParameters.PowerMethod powerMethod = null;
	// quantile of non-centrality distribution - specified if the quantile power
	// method is used
	double quantile = Double.NaN;
	
	// confidence limits for power if requested
	// only available if solving for power in a random design
	ConfidenceInterval confidenceInterval = null;
	
	/**
	 * Create a new GLMMPower object.
	 * 
	 * @param test the statistical test performed
	 * @param alpha the type I error rate
	 * @param nominalPower requested power (for sample size, detectable difference requests)
	 * @param actualPower calculated power
	 * @param sampleSize total sample size
	 * @param betaScale scale factor for beta matrix (roughly interpreted as detectable difference)
	 * @param sigmaScale scale factor for error covariance matrix
	 * @param method method used for power calculation
	 */
	public GLMMPower(GLMMTestFactory.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method)
	{
		this(test, alpha, nominalPower,actualPower,sampleSize,
				betaScale, sigmaScale, method, Double.NaN, null);
	}
	
	/**
	 * Create a new GLMMPower object.
	 * 
	 * @param test the statistical test performed
	 * @param alpha the type I error rate
	 * @param nominalPower requested power (for sample size, detectable difference requests)
	 * @param actualPower calculated power
	 * @param sampleSize total sample size
	 * @param betaScale scale factor for beta matrix (roughly interpreted as detectable difference)
	 * @param sigmaScale scale factor for error covariance matrix
	 * @param method method used for power calculation
	 * @param confidenceInterval confidence interval if requested
	 */
	public GLMMPower(GLMMTestFactory.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method,
			ConfidenceInterval confidenceInterval)
	{
		this(test, alpha, nominalPower,actualPower,sampleSize,
				betaScale, sigmaScale, method, Double.NaN, confidenceInterval);
	}
	
	/**
	 * Create a new GLMMPower object.
	 * 
	 * @param test the statistical test performed
	 * @param alpha the type I error rate
	 * @param nominalPower requested power (for sample size, detectable difference requests)
	 * @param actualPower calculated power
	 * @param sampleSize total sample size
	 * @param betaScale scale factor for beta matrix (roughly interpreted as detectable difference)
	 * @param sigmaScale scale factor for error covariance matrix
	 * @param method method used for power calculation
	 * @param quantile optional quantile value (for quantile power only)
	 */
	public GLMMPower(GLMMTestFactory.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method,
			double quantile)
	{
		this(test, alpha, nominalPower,actualPower,sampleSize,
				betaScale, sigmaScale, method, quantile, null);
	}
	
	/**
	 * Create a new GLMMPower object.
	 * 
	 * @param test the statistical test performed
	 * @param alpha the type I error rate
	 * @param nominalPower requested power (for sample size, detectable difference requests)
	 * @param actualPower calculated power
	 * @param sampleSize total sample size
	 * @param betaScale scale factor for beta matrix (roughly interpreted as detectable difference)
	 * @param sigmaScale scale factor for error covariance matrix
	 * @param method method used for power calculation
	 * @param quantile optional quantile value (for quantile power only)
	 * @param confidenceInterval confidence interval if requested
	 */
	public GLMMPower(GLMMTestFactory.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method,
			double quantile, ConfidenceInterval confidenceInterval)
	{
		super(nominalPower, actualPower, sampleSize, alpha);
		this.test = test;
		this.betaScale = betaScale;
		this.sigmaScale = sigmaScale;
		this.powerMethod = method;
		this.quantile = quantile;
		this.confidenceInterval = confidenceInterval;
	}
	
	/**
	 * Get the regression coefficient scale factor
	 * @return beta scale factor
	 */
	public double getBetaScale()
	{
		return betaScale;
	}

	/**
	 * Set the regression coefficient scale factor
	 * @param betaScale
	 */
	public void setBetaScale(double betaScale)
	{
		this.betaScale = betaScale;
	}

	/**
	 * Get the error covariance scale factor
	 * @return error covariance scale factor
	 */
	public double getSigmaScale()
	{
		return sigmaScale;
	}

	/**
	 * Set error covariance scale factor
	 * 
	 * @param sigmaScale error covariance scale factor
	 */
	public void setSigmaScale(double sigmaScale)
	{
		this.sigmaScale = sigmaScale;
	}
	
	/**
	 * Get the statistical test performed
	 * @return statistical test
	 */
	public GLMMTestFactory.Test getTest()
	{
		return test;
	}

	/**
	 * Set the statistical test
	 * @param test the statistical test
	 * @see edu.cudenver.bios.power.parameters.GLMMPowerParameters.Test
	 */
	public void setTest(GLMMTestFactory.Test test)
	{
		this.test = test;
	}
	
	/**
	 * Get the power calculation method
	 * @return power method
	 * @see edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod 
	 */
	public GLMMPowerParameters.PowerMethod getPowerMethod()
	{
		return powerMethod;
	}

	/**
	 * Set the power calculation method
	 * @param powerMethod
	 * 	@see edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod 
	 */
	public void setPowerMethod(GLMMPowerParameters.PowerMethod powerMethod)
	{
		this.powerMethod = powerMethod;
	}

	/**
	 * Get the quantile used for quantile power.  This will be set to NaN
	 * unless the current power method is quantile power.
	 * @return quantile
	 */
	public double getQuantile()
	{
		return quantile;
	}

	/**
	 * Set the quantile.  This should only be set if the power method is
	 * set to quantile power.
	 * @param quantile the quantile of the non-centrality distribution
	 */
	public void setQuantile(double quantile)
	{
		this.quantile = quantile;
	}

	/**
	 * Get the confidence interval associated with power	
	 * @return confidence interval
	 */
	public ConfidenceInterval getConfidenceInterval()
	{
		return confidenceInterval;
	}

	/**
	 * Set the confidence interval associated with power	
	 * @param confidenceInterval the power confidence interval
	 */
	public void setConfidenceInterval(ConfidenceInterval confidenceInterval)
	{
		this.confidenceInterval = confidenceInterval;
	}

	/**
	 * Output the object as XML
	 */
	public String toXML()
	{
		StringBuffer buffer = new StringBuffer();
		buffer.append("<glmmPower");
		buffer.append(" alpha='");
		buffer.append(alpha);
		buffer.append("' nominalPower='");
		buffer.append(nominalPower);
		buffer.append("' actualPower='");
		buffer.append(actualPower);
		buffer.append("' sampleSize='");
		buffer.append(totalSampleSize);
		buffer.append("' betaScale='");
		buffer.append(betaScale);
		buffer.append("' sigmaScale='");
		buffer.append(sigmaScale);
		buffer.append("' test='");
		buffer.append(test);
		if (powerMethod != null)
		{
			buffer.append("' powerMethod='");
			buffer.append(powerMethod);
		}
		if (quantile != Double.NaN)
		{
			buffer.append("' quantile='");
			buffer.append(quantile);
		}
		buffer.append("' />");
		return buffer.toString();

	}
}
