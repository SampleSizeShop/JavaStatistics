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
	HashMap<String, RealMatrix> twoFactorInteractionContrasts = new HashMap<String, RealMatrix>();
	HashMap<String, RealMatrix> threeFactorInteractionContrasts = new HashMap<String, RealMatrix>();
	
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
	
	public void addTwoFactorInteractionContrast(String factor1Name, String factor2Name, RealMatrix contrast)
	{
		twoFactorInteractionContrasts.put(buildTwoFactorKey(factor1Name, factor2Name), contrast);
	}
	
	public void addThreeFactorInteractionContrast(String factor1Name, String factor2Name, 
			String factor3Name, RealMatrix contrast)
	{
		threeFactorInteractionContrasts.put(buildThreeFactorKey(factor1Name, factor2Name, factor3Name), contrast);
	}
	
	public RealMatrix getMainEffectContrast(String factorName)
	{
		return mainEffectContrasts.get(factorName);
	}
	
	public RealMatrix getTwoFactorInteractionContrast(String factor1Name, String factor2Name)
	{
		return twoFactorInteractionContrasts.get(buildTwoFactorKey(factor1Name, factor2Name));
	}
	
	public RealMatrix getThreeFactorInteractionContrast(String factor1Name, String factor2Name, String factor3Name)
	{
		return threeFactorInteractionContrasts.get(buildThreeFactorKey(factor1Name, factor2Name, factor3Name));
	}
	
	/**
	 * Combine the two factor names into a single hash map key
	 * @param factor1Name name of first factor
	 * @param factor2Name name of second factor
	 * @return hashmap key for two factor interaction
	 */
	private String buildTwoFactorKey(String factor1Name, String factor2Name)
	{
		StringBuffer buffer = new StringBuffer();
		buffer.append("F1=");
		buffer.append(factor1Name);
		buffer.append("&F2=");
		buffer.append(factor2Name);
		return buffer.toString();
	}
	
	/**
	 * Combine the three factor names into a single hash map key
	 * @param factor1Name name of first factor
	 * @param factor2Name name of second factor
	 * @return hashmap key for three factor interaction
	 */
	private String buildThreeFactorKey(String factor1Name, String factor2Name, String factor3Name)
	{
		StringBuffer buffer = new StringBuffer();
		buffer.append("F1=");
		buffer.append(factor1Name);
		buffer.append("&F2=");
		buffer.append(factor2Name);
		buffer.append("&F3=");
		buffer.append(factor3Name);
		return buffer.toString();
	}
}
