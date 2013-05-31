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
import java.util.List;

import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.OrthogonalPolynomialContrast.ContrastType;
import edu.cudenver.bios.utils.Factor;

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
	ArrayList<OrthogonalPolynomialContrast> contrastList = new ArrayList<OrthogonalPolynomialContrast>();
	
	/**
	 * Create an empty contrast collection
	 */
	public OrthogonalPolynomialContrastCollection() {}
	
	/**
	 * Add a contrast to the collection
	 * @param contrast the contrast
	 */
	public void addContrast(OrthogonalPolynomialContrast contrast)
	throws IllegalArgumentException
	{
		if (findContrast(contrast) != null)
			throw new IllegalArgumentException("Contrast already in collection");
		contrastList.add(contrast);
		int foo = 1;
	}
	
	/**
	 * Get the contrast to test the grand mean
	 * @return grand mean contrast
	 */
	public OrthogonalPolynomialContrast getGrandMean()
	{
		for(OrthogonalPolynomialContrast contrast: contrastList)
		{
			if (ContrastType.GRAND_MEAN == contrast.getType())
				return contrast;
		}
		return null;
	}
	
	/**
	 * Check if a contrast already exists in the collection for the specified
	 * type and list of factors
	 * 
	 * @param type contrast type
	 * @param factorList list of factors for the contrast (null if grand mean)
	 * @return if found, the matching contrast.  Null if not found.
	 */
	private OrthogonalPolynomialContrast findContrast(OrthogonalPolynomialContrast compContrast)
	{
		for(OrthogonalPolynomialContrast contrast: contrastList)
		{
			if (contrast.compareTo(compContrast) == 0) return contrast;
		}
		return null;
	}
	
	public OrthogonalPolynomialContrast getMainEffectContrast(Factor factor)
	{
		ArrayList<Factor> factorList = new ArrayList<Factor>(1);
		factorList.add(factor);
		OrthogonalPolynomialContrast compContrast = 
			new OrthogonalPolynomialContrast(ContrastType.MAIN_EFFECT, factorList, null);
		for(OrthogonalPolynomialContrast contrast: contrastList)
		{
			if (contrast.compareTo(compContrast) == 0) return contrast;
		}
		return null;
	}
	
	public OrthogonalPolynomialContrast getInteractionContrast(List<Factor> factorList)
	{
		OrthogonalPolynomialContrast compContrast = 
			new OrthogonalPolynomialContrast(ContrastType.INTERACTION, factorList, null);
		for(OrthogonalPolynomialContrast contrast: contrastList)
		{
			if (contrast.compareTo(compContrast) == 0) return contrast;
		}
		return null;
	}
	
	public List<OrthogonalPolynomialContrast> getContrastList()
	{
		return contrastList;
	}
}
