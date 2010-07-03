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

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

/**
 * Pojo containing a description of the power result
 */
public class GLMMPower extends Power
{
	// scale factor for beta matrix
	double betaScale;
	// scale factor for the sigma error matrix
	double sigmaScale;
	// statistical test performed
	GLMMPowerParameters.Test test;
	
	/* optional parameters */
	// power method - specified if a baseline covariate is used
	GLMMPowerParameters.PowerMethod powerMethod = null;
	// quantile of non-centrality distribution - specified if the quantile power
	// method is used
	double quantile = Double.NaN;
	
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
	public GLMMPower(GLMMPowerParameters.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method)
	{
		super(nominalPower, actualPower, sampleSize, alpha);
		this.test = test;
		this.betaScale = betaScale;
		this.sigmaScale = sigmaScale;
		this.powerMethod = method;
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
	public GLMMPower(GLMMPowerParameters.Test test, double alpha, 
			double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale, 
			GLMMPowerParameters.PowerMethod method,
			double quantile)
	{
		super(nominalPower, actualPower, sampleSize, alpha);
		this.test = test;
		this.betaScale = betaScale;
		this.sigmaScale = sigmaScale;
		this.powerMethod = method;
		this.quantile = quantile;
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
	public GLMMPowerParameters.Test getTest()
	{
		return test;
	}

	/**
	 * Set the statistical test
	 * @param test the statistical test
	 * @see GLMMPowerParameters.Test
	 */
	public void setTest(GLMMPowerParameters.Test test)
	{
		this.test = test;
	}
	
	/**
	 * Get the power calculation method
	 * @return power method
	 * @see GLMMPowerParameters.PowerMethod 
	 */
	public GLMMPowerParameters.PowerMethod getPowerMethod()
	{
		return powerMethod;
	}

	/**
	 * Set the power calculation method
	 * @param powerMethod
	 * 	@see GLMMPowerParameters.PowerMethod 
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
