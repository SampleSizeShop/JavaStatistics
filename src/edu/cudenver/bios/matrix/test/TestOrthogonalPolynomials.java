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
package edu.cudenver.bios.matrix.test;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.OrthogonalPolynomialContrastCollection;
import edu.cudenver.bios.matrix.OrthogonalPolynomials;
import edu.cudenver.bios.utils.Factor;
import junit.framework.TestCase;

/**
 * Test cases for OrthogonalPolynomials class
 * @author Sarah Kreidler
 *
 */
public class TestOrthogonalPolynomials extends TestCase
{	
	/**
	 * Test a one-factor U matrix 
	 */
	public void testOneFactorContrast()
	{
		double[] x = {1,2,3,4};
		String name = "x";
		
		ArrayList<Factor> factorList = new ArrayList<Factor>();
		factorList.add(new Factor(name, x));
		OrthogonalPolynomialContrastCollection contrastCollection = 
			OrthogonalPolynomials.withinSubjectContrast(factorList);
	
		printMatrix("Grand mean", contrastCollection.getGrandMean());
		printMatrix("One factor contrast", contrastCollection.getMainEffectContrast(name));		
	}

	/**
	 * Test a two-factor U matrix 
	 */
	public void testTwoFactorContrast()
	{
		double[] x1 = {1,2,4};
		String x1Name = "x1";
		double[] x2 = {1,3,5};
		String x2Name = "x2";
		
		ArrayList<Factor> factorList = new ArrayList<Factor>();
		factorList.add(new Factor(x1Name, x1));
		factorList.add(new Factor(x2Name, x2));
		OrthogonalPolynomialContrastCollection contrastCollection = 
			OrthogonalPolynomials.withinSubjectContrast(factorList);
	
		printMatrix("Grand mean", contrastCollection.getGrandMean());
		printMatrix("Main effect x1", contrastCollection.getMainEffectContrast(x1Name));
		printMatrix("Main effect x2", contrastCollection.getMainEffectContrast(x2Name));
		
		ArrayList<String> names = new ArrayList<String>(2);
		names.add(x1Name);
		names.add(x2Name);
		printMatrix("Two factor contrast", contrastCollection.getInteractionContrast(names));		
	}
	
	/**
	 * Test a three-factor U matrix 
	 */
	public void testThreeFactorContrast()
	{
		double[] x1 = {1,2,3};
		String x1Name = "x1";
		double[] x2 = {1,2,3};
		String x2Name = "x2";
		double[] x3 = {1,2,3};
		String x3Name = "x3";
		ArrayList<Factor> factorList = new ArrayList<Factor>();
		factorList.add(new Factor(x1Name, x1));
		factorList.add(new Factor(x2Name, x2));
		factorList.add(new Factor(x3Name, x3));
		OrthogonalPolynomialContrastCollection contrastCollection = 
			OrthogonalPolynomials.withinSubjectContrast(factorList);
	
		printMatrix("Grand mean", contrastCollection.getGrandMean());
		printMatrix("Main effect x1", contrastCollection.getMainEffectContrast(x1Name));
		printMatrix("Main effect x2", contrastCollection.getMainEffectContrast(x2Name));
		printMatrix("Main effect x3", contrastCollection.getMainEffectContrast(x3Name));
		
		ArrayList<String> names = new ArrayList<String>(5);
		names.add(x1Name);
		names.add(x2Name);
		printMatrix("Two factor contrast: x1, x2", contrastCollection.getInteractionContrast(names));		
		names.clear();
		names.add(x1Name);
		names.add(x3Name);
		printMatrix("Two factor contrast: x1, x3", contrastCollection.getInteractionContrast(names));		
		names.clear();
		names.add(x2Name);
		names.add(x3Name);
		printMatrix("Two factor contrast: x2, x3", contrastCollection.getInteractionContrast(names));		
		names.clear();
		names.add(x1Name);
		names.add(x2Name);
		names.add(x3Name);
		printMatrix("Three factor contrast", contrastCollection.getInteractionContrast(names));		
	}
	
	/**
	 * Write the matrix to std out
	 * @param m
	 */
	private void printMatrix(String title, RealMatrix m)
	{
		System.out.println(title);
	    DecimalFormat Number = new DecimalFormat("#0.000");
		for(int row = 0; row < m.getRowDimension(); row++)
		{
			for(int col= 0; col < m.getColumnDimension(); col++)
			{
				System.out.print(Number.format(m.getEntry(row, col)) + "\t");
			}
			System.out.print("\n");
		}
	}
}
