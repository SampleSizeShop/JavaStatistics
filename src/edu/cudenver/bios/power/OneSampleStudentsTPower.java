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

/**
 * Power results for the One Sample Student's T Test
 * @author Sarah Kreidler
 *
 */
public class OneSampleStudentsTPower extends Power
{
	double sigma;
	double detectableDifference;
	
	/**
	 * Constructor
	 * @param alpha type I error
	 * @param nominalPower target power
	 * @param actualPower calculated power
	 * @param sampleSize total sample size
	 * @param detectableDifference effect size
	 * @param sigma
	 */
	public OneSampleStudentsTPower(double alpha, double nominalPower, double actualPower,
			int sampleSize, double detectableDifference, double sigma)
	{
		super(nominalPower, actualPower, sampleSize, alpha);
		this.sigma = sigma;
		this.detectableDifference = detectableDifference;
	}
	
	/**
	 * Get variance used in power calculation
	 * @return sigma
	 */
	public double getSigma()
	{
		return sigma;
	}

	/**
	 * Set the variance 
	 * @param sigma
	 */
	public void setSigma(double sigma)
	{
		this.sigma = sigma;
	}

	/**
	 * Get the detectable difference
	 * @return detectable difference
	 */
	public double getDetectableDifference()
	{
		return detectableDifference;
	}

	/**
	 * Set the detectable difference
	 * @param detectableDifference
	 */
	public void setDetectableDifference(double detectableDifference)
	{
		this.detectableDifference = detectableDifference;
	}

	/**
	 * Return power result as XML
	 */
	public String toXML()
	{
		return "<power alpha='"+ alpha +"' nominalPower='" + nominalPower + "' actualPower='"+ actualPower +
			"' sampleSize='"+ totalSampleSize + "' difference='"+detectableDifference +"' sigma='"+ sigma +"' />";
	}
}
