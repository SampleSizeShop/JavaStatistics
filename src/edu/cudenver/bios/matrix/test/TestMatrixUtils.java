package edu.cudenver.bios.matrix.test;

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
import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.MatrixUtils;

import junit.framework.TestCase;

public class TestMatrixUtils extends TestCase
{
	private static final double TOLERANCE = 1.0E-15;
	double[][] dataA = {{1,2},{3,4}};
	double[][] dataB = {{5,6},{7,8}};
	double[][] dataC = {{9,10},{11,12}};
	double[][] dataD = {{13,14,15},{16,17,18}};
	Array2DRowRealMatrix matrixA = new Array2DRowRealMatrix(dataA);
	Array2DRowRealMatrix matrixB = new Array2DRowRealMatrix(dataB);
	Array2DRowRealMatrix matrixC = new Array2DRowRealMatrix(dataC);
	Array2DRowRealMatrix matrixD = new Array2DRowRealMatrix(dataD);
	
	public void testKroneckerTwoMatrix()
	{
		double[][] data = {{5,6,10,12},{7,8,14,16},{15,18,20,24},{21,24,28,32}};
		RealMatrix correct = new Array2DRowRealMatrix(data);
		RealMatrix result = MatrixUtils.getKroneckerProduct(matrixA, matrixB);
		assertTrue(compareMatrices(correct, result));

	}
	
	public void testKroneckerTwoMatrixListEntry()
	{
		double[][] data = {{5,6,10,12},{7,8,14,16},{15,18,20,24},{21,24,28,32}};
		RealMatrix correct = new Array2DRowRealMatrix(data);
		
		ArrayList<RealMatrix> list = new ArrayList<RealMatrix>();
		list.add(matrixA);
		list.add(matrixB);
		RealMatrix result = MatrixUtils.getKroneckerProduct(list);
		
		assertTrue(compareMatrices(correct, result));

	}
	
	public void testKroneckerThreeMatrixEqualDimension()
	{
		double[][] data = {
				{45,50,54,60,90,100,108,120},
				{55,60,66,72,110,120,132,144},
				{63,70,72,80,126,140,144,160},
				{77,84,88,96,154,168,176,192},
				{135,150,162,180,180,200,216,240},
				{165,180,198,216,220,240,264,288},
				{189,210,216,240,252,280,288,320},
				{231,252,264,288,308,336,352,384}
				};
		RealMatrix correct = new Array2DRowRealMatrix(data);
		
		ArrayList<RealMatrix> list = new ArrayList<RealMatrix>();
		list.add(matrixA);
		list.add(matrixB);
		list.add(matrixC);
		RealMatrix result = MatrixUtils.getKroneckerProduct(list);
		
		assertTrue(compareMatrices(correct, result));
	}
	
	public void testKroneckerThreeMatrixUnequalDimension()
	{
		double[][] data = {
				{65,70,75,78,84,90,130,140,150,156,168,180},
				{80,85,90,96,102,108,160,170,180,192,204,216},
				{91,98,105,104,112,120,182,196,210,208,224,240},
				{112,119,126,128,136,144,224,238,252,256,272,288},
				{195,210,225,234,252,270,260,280,300,312,336,360},
				{240,255,270,288,306,324,320,340,360,384,408,432},
				{273,294,315,312,336,360,364,392,420,416,448,480},
				{336,357,378,384,408,432,448,476,504,512,544,576}
				};
		RealMatrix correct = new Array2DRowRealMatrix(data);
		
		ArrayList<RealMatrix> list = new ArrayList<RealMatrix>();
		list.add(matrixA);
		list.add(matrixB);
		list.add(matrixD);
		RealMatrix result = MatrixUtils.getKroneckerProduct(list);
		
		assertTrue(compareMatrices(correct, result));
	}
	
	/**
	 * Determines if dimensions and all cells are equal between the
	 * two matrices
	 * 
	 * @param A matrix A
	 * @param B matrix B
	 * @return true if equal, false if not equal
	 */
	private boolean compareMatrices(RealMatrix A, RealMatrix B)
	{
		if (A.getRowDimension() == B.getRowDimension() &&
				A.getColumnDimension() == B.getColumnDimension())
		{
			for(int r = 0; r < A.getRowDimension(); r++)
			{
				for(int c = 0; c < A.getColumnDimension(); c++)
				{
					if (Double.compare(A.getEntry(r, c), B.getEntry(r, c)) != 0) 
					{
						System.err.println("Value mismatch at row=" + r + ", column=" + c);
						return false;
					}
				}
			}
		}
		else
		{
			System.err.println("Dimensions do not match: A ("+ 
					A.getRowDimension() + " x " + A.getColumnDimension() + 
					"), B ("+ B.getRowDimension() +" x "+ 
					B.getColumnDimension()+ ")");
			return false;
		}
		return true;
	}
		
}
