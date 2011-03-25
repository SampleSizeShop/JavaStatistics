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

import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.OrthogonalPolynomials;
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
		RealMatrix contrast = OrthogonalPolynomials.withinSubjectContrast(x);
	
		printMatrix("One factor contrast", contrast);		
	}

	/**
	 * Test a two-factor U matrix 
	 */
	public void testTwoFactorContrast()
	{
		double[] x1 = {1,2,4};
		double[] x2 = {1,3,5};
		RealMatrix contrast = OrthogonalPolynomials.withinSubjectContrast(x1,x2);
	
		printMatrix("Two factor contrast", contrast);		
	}
	
	/**
	 * Test a three-factor U matrix 
	 */
	public void testThreeFactorContrast()
	{
		double[] x1 = {1,2,3};
		double[] x2 = {1,2,3};
		double[] x3 = {1,2,3};
		RealMatrix contrast = OrthogonalPolynomials.withinSubjectContrast(x1,x2,x3);
	
		printMatrix("Three factor contrast", contrast);		
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
