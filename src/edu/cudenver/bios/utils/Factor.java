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
 * Container class for factor information, which includes a label
 * and a list of values for the factor.  For repeated measures, the
 * values would represent the spacing of measurements in time. 
 * For cluster samples, the values are likely sequential assuming
 * exchangeable samples (for example, a cluster of 5 schools would
 * have values 1,2,3,4,5).
 * 
 * @author Sarah Kreidler
 *
 */
public class Factor
{
	// factor name
	String name;
	// spacing of repeated measures
	double[] values;
	
	/**
	 * Constructor
	 * @param name name of the factor
	 * @param values spacing of repeated measures
	 */
	public Factor(String name, double[] values)
	{
		this.name = name;
		this.values = values;
	}

	/**
	 * Get the name of the factor
	 * @return factor name
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * Get the factor values
	 * @return factor value array
	 */
	public double[] getValues()
	{
		return values;
	}
	
	
}
