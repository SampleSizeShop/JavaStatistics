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
package edu.cudenver.bios.matrix;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math.linear.RealMatrix;

/**
 * A collection of orthogonal polynomial contrasts for up to 3 factors
 * which includes main effects for each factor, pair-wise interactions,
 * and 3-way interaction contrasts.
 * 
 * @author Sarah Kreidler
 *
 */
public class OrthogonalPolynomialContrastCollection
{
	RealMatrix grandMean;
	HashMap<String, RealMatrix> mainEffectContrasts = new HashMap<String, RealMatrix>();
	HashMap<String, RealMatrix> interactionContrasts = new HashMap<String, RealMatrix>();
	
	/**
	 * Create an empty contrast collection
	 */
	public OrthogonalPolynomialContrastCollection() {}
	
	public void setGrandMean(RealMatrix grandMean)
	{
		this.grandMean = grandMean;
	}
	
	public RealMatrix getGrandMean()
	{
		return this.grandMean;
	}
	
	public void addMainEffectContrast(String factorName, RealMatrix contrast)
	{
		mainEffectContrasts.put(factorName, contrast);
	}
	
	public void addInteractionContrast(ArrayList<String> factorNames, RealMatrix contrast)
	{
		interactionContrasts.put(buildFactorKey(factorNames), contrast);
	}
	
	public RealMatrix getMainEffectContrast(String factorName)
	{
		return mainEffectContrasts.get(factorName);
	}
	
	public RealMatrix getInteractionContrast(ArrayList<String> factorNames)
	{
		return interactionContrasts.get(buildFactorKey(factorNames));
	}
	
	/**
	 * Combine the two factor names into a single hash map key
	 * @param factorNames list of factor names in order
	 * @return hashmap key for interaction contrast
	 */
	private String buildFactorKey(ArrayList<String> factorNames)
	{
		StringBuffer buffer = new StringBuffer();
		int count = 0;
		for(String name : factorNames)
		{
			buffer.append("&F");
			buffer.append(count);
			buffer.append("=");
			buffer.append(name);
			count++;
		}
		return buffer.toString();
	}

}
