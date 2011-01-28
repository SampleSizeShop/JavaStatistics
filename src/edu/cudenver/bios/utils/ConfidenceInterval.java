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
package edu.cudenver.bios.utils;

/**
 * POJO describing a confidence interval
 * 
 * @author Sarah Kreidler
 *
 */
public class ConfidenceInterval
{
	protected double alphaLower;
	protected double alphaUpper;
	protected double lowerLimit;
	protected double upperLimit;
	
	public ConfidenceInterval(double lowerLimit, double upperLimit, double alphaLower,
			double alphaUpper)
	throws IllegalArgumentException
	{
		if (lowerLimit > upperLimit)
			throw new IllegalArgumentException("invalid bounds for confidence interval");
		if (alphaLower <= 0 || alphaLower >= 1 || 
				alphaUpper <= 0 || alphaUpper >= 1 ||
				alphaUpper + alphaLower >= 1)
			throw new IllegalArgumentException("invalid alpha values for confidence interval");
		
		this.alphaLower = alphaLower;
		this.alphaUpper = alphaUpper;
		this.lowerLimit = lowerLimit;
		this.upperLimit = upperLimit;
		
	}
	
	public ConfidenceInterval() {}

	public double getAlphaLower()
	{
		return alphaLower;
	}

	public void setAlphaLower(double alphaLower)
	{
		this.alphaLower = alphaLower;
	}

	public double getAlphaUpper()
	{
		return alphaUpper;
	}

	public void setAlphaUpper(double alphaUpper)
	{
		this.alphaUpper = alphaUpper;
	}

	public double getLowerLimit()
	{
		return lowerLimit;
	}

	public void setLowerLimit(double lowerLimit)
	{
		this.lowerLimit = lowerLimit;
	}

	public double getUpperLimit()
	{
		return upperLimit;
	}

	public void setUpperLimit(double upperLimit)
	{
		this.upperLimit = upperLimit;
	}
	
	public double getConfidenceCoefficient()
	{
		return 1 - alphaUpper - alphaLower;
	}
	
}
