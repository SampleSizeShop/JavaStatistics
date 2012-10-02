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
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.utils.Factor;

public class OrthogonalPolynomialContrast implements Comparable<OrthogonalPolynomialContrast>
{
	// the type of contrast
	public enum ContrastType
	{
		GRAND_MEAN,
		MAIN_EFFECT,
		INTERACTION
	};
	// the factors associated with the contrast
	protected ArrayList<Factor> factorList = new ArrayList<Factor>();
	// the contrast matrix
	protected RealMatrix contrastMatrix = null;
	// the type of contrast
	protected ContrastType type;
	
	/**
	 * Create an orthogonal polynomial contrast object.  A polynomial contrast can be used to
	 * test the grand mean, main-effects or interactions.
	 * 
	 * @param type the type of hypothesis being tested by the 
	 * contrast (grand mean, main-effect, interaction)
	 * @param factorList list of factors tested by the contrast
	 * @param contrastMatrix matrix of contrast coefficients
	 */
	public OrthogonalPolynomialContrast(ContrastType type, List<Factor> factorList, RealMatrix contrastMatrix)
	{
		this.type = type;
		for(Factor factor: factorList) this.factorList.add(factor);
		Collections.sort(this.factorList);
		this.contrastMatrix = contrastMatrix;
	}
	
	/**
	 * Create a grand mean contrast using orthogonal polynomials.
	 * 
	 * @param type the type of hypothesis being tested by the 
	 * contrast (grand mean, main-effect, interaction)
	 * @param factorList list of factors tested by the contrast
	 * @param contrastMatrix matrix of contrast coefficients
	 */
	public OrthogonalPolynomialContrast(RealMatrix contrastMatrix)
	{
		this.type = ContrastType.GRAND_MEAN;
		this.factorList = null;
		this.contrastMatrix = contrastMatrix;
	}
	
	/**
	 * Get the contrast type: grand mean, main effect, or interaction
	 * @return contrast type
	 */
	public ContrastType getType()
	{
		return type;
	}
	
	/**
	 * Get the list of factors for the contrast
	 * @return list of factors for the current contrast.  Null for grand mean contrasts
	 */
	public List<Factor> getFactorList()
	{
		return factorList;
	}
	
	/**
	 * Get the contrast matrix
	 * @return contrast matrix
	 */
	public RealMatrix getContrastMatrix()
	{
		return contrastMatrix;
	}

	@Override
	public int compareTo(OrthogonalPolynomialContrast compContrast)
	{
		// compare the types - top level sort in order: grand mean, main effect, interaction
		if (type != compContrast.getType())
		{
			switch (type)
			{
			case GRAND_MEAN:
				return -1;
			case MAIN_EFFECT:
				if (compContrast.getType() == ContrastType.GRAND_MEAN)
					return 1;
				else
					return -1;
			default:
				return 1;
			}
		}
		else
		{
			// compare factor lists - assumes both are sorted
			List<Factor> compFactorList = compContrast.getFactorList();
			if (factorList.size() != compFactorList.size())
			{
				// sort contrasts with fewer factors first
				return factorList.size() - compFactorList.size();
			}
			else
			{
				// #factors is equal, assume both are sorted and compare
				for(int i = 0; i < factorList.size(); i++)
				{
					int comp = factorList.get(i).compareTo(compFactorList.get(i));
					if (comp != 0) return comp;
				}
			}
		}
		return 0;
	}
}
