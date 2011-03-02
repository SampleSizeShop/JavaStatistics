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

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

public class MatrixUtils
{
	/**
	 * Creates a matrix of equal dimension but with all non-diagonal 
	 * elements set to 0
	 * 
	 * @param matrix the matrix from which to extract the diagonal
	 * @return diagonal matrix
	 */
	public static RealMatrix toDiagonalMatrix(RealMatrix matrix)
	throws IllegalArgumentException
	{
		if (matrix == null)
			throw new IllegalArgumentException("null input matrix");
		
		double[][] zData = new double[matrix.getRowDimension()][matrix.getColumnDimension()]; 
		for(int row = 0; row < matrix.getRowDimension(); row++)
		{
			for(int col = 0; col < matrix.getColumnDimension(); col++)
			{
				if (row == col) 
					zData[row][col] = matrix.getEntry(row, col);
				else
					zData[row][col] = 0;
			}
		}
		return new Array2DRowRealMatrix(zData);
	}
	
	public static RealMatrix horizontalDirectProduct(RealMatrix matrix1, RealMatrix matrix2)
	throws IllegalArgumentException
	{
		if (matrix1 == null ||matrix2 == null)
			throw new IllegalArgumentException("null input matrix");
		if (matrix1.getRowDimension() != matrix2.getRowDimension())
			throw new IllegalArgumentException("input matrices must have equal row dimension");
		
		int mRows = matrix1.getRowDimension();
		int m1Cols = matrix1.getColumnDimension();
		int m2Cols = matrix2.getColumnDimension();
		
		double[][] productData = new double[mRows][m2Cols * m2Cols]; 
		RealMatrix productMatrix = new Array2DRowRealMatrix(productData);
		for(int col = 0; col < m1Cols; col++)
		{
			for(int row = 0; row < mRows; row++)
			{
				RealMatrix m2Row = matrix2.getRowMatrix(row);
				productMatrix.setSubMatrix((m2Row.scalarMultiply(matrix1.getEntry(row, col))).getData(), 
						row, col * m2Cols);
			}
		}
		return productMatrix;
	}
	
	public static RealMatrix KroneckerProduct(RealMatrix matrix1, RealMatrix matrix2)
	{
		if (matrix1 == null ||matrix2 == null)
			throw new IllegalArgumentException("null input matrix");
		
		int m1Rows = matrix1.getRowDimension();
		int m1Cols = matrix1.getColumnDimension();
		int m2Rows = matrix2.getRowDimension();
		int m2Cols = matrix2.getColumnDimension();
		
		double[][] productData = new double[m1Rows*m2Rows][m1Cols * m2Cols]; 
		RealMatrix productMatrix = new Array2DRowRealMatrix(productData);
		for(int col = 0; col < m1Cols; col++)
		{
			for(int row = 0; row < m1Rows; row++)
			{
				productMatrix.setSubMatrix((matrix2.scalarMultiply(matrix1.getEntry(row, col))).getData(), 
						row * m2Rows, col * m2Cols);
			}
		}
		
		return productMatrix;
	}
}
