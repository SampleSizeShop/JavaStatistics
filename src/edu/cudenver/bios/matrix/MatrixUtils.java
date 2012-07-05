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

import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import jsc.distributions.Normal;

public class MatrixUtils
{
	
	/**
	 * Creates a matrix of equal dimension but with all non-diagonal 
	 * elements set to 0
	 * 
	 * @param matrix the matrix from which to extract the diagonal
	 * @return diagonal matrix
	 */
	public static RealMatrix getDiagonalMatrix(RealMatrix matrix)
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
	
	/**
	 * Calculate the horizontal direct product of two matrices
	 * 
	 * @param matrix1 first matrix
	 * @param matrix2 second matrix
	 * @return horizontal direct product of matrix 1 and matrix 2
	 */
	public static RealMatrix getHorizontalDirectProduct(RealMatrix matrix1, RealMatrix matrix2)
	throws IllegalArgumentException
	{
		if (matrix1 == null || matrix2 == null)
			throw new IllegalArgumentException("null input matrix");
		if (matrix1.getRowDimension() != matrix2.getRowDimension())
			throw new IllegalArgumentException("input matrices must have equal row dimension");
		
		int mRows = matrix1.getRowDimension();
		int m1Cols = matrix1.getColumnDimension();
		int m2Cols = matrix2.getColumnDimension();
		
		double[][] productData = new double[mRows][m1Cols * m2Cols]; 
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
	
	/**
	 * Calculate the Kronecker product of a list of matrices.  Applies the
	 * Kronecker product following the list order (i.e. from left to right).
	 * 
	 * @param matrixList list of matrices
	 * @return Kronecker product of matrices
	 */
	public static RealMatrix getKroneckerProduct(List<RealMatrix> matrixList)
	{
		if (matrixList == null || matrixList.size() <= 0)
			throw new IllegalArgumentException("no input matrices");
		if (matrixList.size() == 1)
			return matrixList.get(0); // nothing to do, only one matrix
		
		// calculate the dimensions of the Kronecker product matrix
		int totalRows = 1;
		int totalCols = 1;
		for(RealMatrix matrix: matrixList)
		{
			totalRows *= matrix.getRowDimension();
			totalCols *= matrix.getColumnDimension();
		}
		
		// create a matrix to hold the data
		double[][] productData = new double[totalRows][totalCols]; 
		// initialize to 1 (allows us to multiple the contents of each matrix 
		// onto the result sequentially
		for(int prodRow = 0; prodRow < totalRows; prodRow++)
		{
			for(int prodCol = 0; prodCol < totalCols; prodCol++)
			{
				productData[prodRow][prodCol] = 1;
			}
		}
		
		// multiply the contents of each matrix onto the result
		int maxRow = totalRows;
		int maxCol = totalCols;
		for(RealMatrix matrix: matrixList)
		{
			maxRow /= matrix.getRowDimension();
			maxCol /= matrix.getColumnDimension();
			int matrixRow = 0;
			int matrixCol = 0;
			// multiply onto the result
			for(int prodRow = 0, sectionRow = 0; prodRow < totalRows; prodRow++, sectionRow++)
			{
				matrixCol = 0;
				double value = matrix.getEntry(matrixRow, matrixCol);
				for(int prodCol = 0, sectionCol = 0; prodCol < totalCols; prodCol++, sectionCol++)
				{
					productData[prodRow][prodCol] *= value;
					if (sectionCol >= maxCol - 1) 
					{
						matrixCol++;
						if (matrixCol >= matrix.getColumnDimension()) matrixCol = 0;
						sectionCol = -1;
						value = matrix.getEntry(matrixRow, matrixCol);
					}
				}
				if (sectionRow >= maxRow-1)
				{
					matrixRow++;
					if (matrixRow >= matrix.getRowDimension()) matrixRow = 0;
					sectionRow = -1;
				}
			}	
		}
		
		// return a new matrix containing the Kronecker product
		return new Array2DRowRealMatrix(productData);
	}
	
	/**
	 * Calculate the Kronecker product of two matrices
	 * 
	 * @param matrix1 first matrix
	 * @param matrix2 second matrix
	 * @return Kronecker product of matrix 1 and matrix 2
	 */
	public static RealMatrix getKroneckerProduct(RealMatrix matrix1, RealMatrix matrix2)
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
	
	/**
	 * Creates a matrix of the specified size with all cells set to the specified value
	 * 
	 * @param rows row dimension
	 * @param cols column dimension
	 * @param value fill value
	 * @return rows x cols matrix with all cells set to the specified value
	 */
	public static RealMatrix getRealMatrixWithFilledValue(int rows, int cols, double value)
	{
		double[][] data = new double[rows][cols];
		
		for(int r = 0; r < rows; r++)
		{
			for(int c = 0; c < cols; c++)
			{
				data[r][c] = value;
			}
		}
		return new Array2DRowRealMatrix(data);
	}
	
	/**
	 * Creates a matrix of the specified size with all cells set to the specified value
	 * 
	 * @param rows row dimension
	 * @param cols column dimension
	 * @param value fill value
	 * @return rows x cols matrix with all cells set to the specified value
	 */
	public static RealMatrix getRandomColumnMatrix(int rows, double mean, double variance, int seed)
	{
		double[] data = new double[rows];
        Normal dist = new Normal(mean, Math.sqrt(variance));
        dist.setSeed(seed);
        
        for(int r = 0; r < rows; r++)
        {
        	data[r] = dist.random();
        }
		return new Array2DRowRealMatrix(data);
	}
	
	/**
	 * Horizontally append two matrices
	 * @param matrix
	 * @param column
	 * @return the combined matrix
	 * @throws IllegalArgumentException
	 */
	public static RealMatrix getHorizontalAppend(RealMatrix m1, RealMatrix m2)
	throws IllegalArgumentException
	{
		if (m1 == null || m2 == null)
			throw new IllegalArgumentException("Missing required argument");
		if (m1.getRowDimension() != m2.getRowDimension())
			throw new IllegalArgumentException("Row dimensions must be equal");
		
		RealMatrix newMatrix = 
			new Array2DRowRealMatrix(m1.getRowDimension(),
					m1.getColumnDimension()+m2.getColumnDimension());
		newMatrix.setSubMatrix(m1.getData(), 0, 0);
		newMatrix.setSubMatrix(m2.getData(), 0, m1.getColumnDimension());
		
		return newMatrix;		
	}
	
	/**
	 * Regenerates the design matrix assuming a fixed (GLMM(F)) design
	 * 
	 * @return full design matrix
	 */
	public static RealMatrix getFullDesignMatrix(RealMatrix designEssence, int sampleSize)
	{
		if( designEssence == null){
			throw new IllegalArgumentException("designEssence is null.");
		}
		return MatrixUtils.getKroneckerProduct(designEssence, 
				MatrixUtils.getRealMatrixWithFilledValue(sampleSize, 1, 1));
	}
	
	/**
	 * Regenerates the design matrix and fills any random columns with a new
	 * realization of random values based on a normal distribution.
	 * 
	 * @return full design matrix
	 */
	public static RealMatrix getFullDesignMatrix(RealMatrix designEssence, 
			RealMatrix sigmaGaussianRandom, int sampleSize, int seed)
	{
		RealMatrix fixed = MatrixUtils.getKroneckerProduct(designEssence, 
				MatrixUtils.getRealMatrixWithFilledValue(sampleSize, 1, 1));
			// GLMM(F, g), so we need to append a random column
			return MatrixUtils.getHorizontalAppend(fixed, 
					MatrixUtils.getRandomColumnMatrix(fixed.getRowDimension(), 0, 
							sigmaGaussianRandom.getEntry(0, 0), seed));
	}
	
    /**
     * Returns the total sample size for the current per group sample size
     */
    public static int getTotalSampleSize(RealMatrix designEssence, int perGroupSize)
    {
    	return designEssence.getRowDimension() * perGroupSize;
    }
    
    /**
     * Return the element wise sum of squares
     * @param matrix input matrix
     * @return sum of squares
     */
    public static double getSumOfSquares(RealMatrix matrix)
    {
    	if( matrix == null){
    		throw new IllegalArgumentException("Null matrix not allowed for " +
    				"getSumOfSquares()");
    	}
    	double sum = 0.0;
		for(int r = 0; r < matrix.getRowDimension(); r++)
		{
			for(int c = 0; c < matrix.getColumnDimension(); c++)
			{
				double value = matrix.getEntry(r, c);
				sum += value * value;
			}
		}
		return sum;
    }
    
    /**
     * Returns the vec(matrix).
     * @param RealMatrix matrix 
     * @return RealMatrix representing vec(matrix).
     */
    public static RealMatrix getVecMatrix(RealMatrix matrix){
    	if( matrix == null){
    		throw new IllegalArgumentException("Null matrix not allowed for " +
    				"getVecMatrix()");
    	}
    	int numRows = matrix.getRowDimension();
    	int numCols = matrix.getColumnDimension();
    	RealMatrix vec = new Array2DRowRealMatrix( numRows*numCols, 1 );
    	int newRowNum = 0;
    	
    	//loop through each column
    	for(int c = 0; c < numCols; c++){
    	    //insert column values into new r x 1 matrix
    		for( int r = 0; r < numRows; r++, newRowNum++){
    			vec.setEntry(newRowNum, 0, matrix.getEntry(r, c) );
    	    }
    	}
    	return vec;
    }
    
    /**
     * This method will return the element-wise product of matrixA
     * and matrixB.
     * In order to perform element-wise multiplication, the incoming
     * matrices must have the same dimensions.
     * @param matrixA
     * @param matrixB
     * @return RealMatrix which is the element-wise product of A and B.
     */
    public static RealMatrix getElementWiseProduct(RealMatrix matrixA, RealMatrix matrixB)
    throws IllegalArgumentException{
    	if( matrixA == null ||
    		matrixB == null ||
    		! areDimensionsEqual(matrixA, matrixB) ){
    		throw new IllegalArgumentException("Both matrices must be non-null " +
    				"and matrix dimensions must be equal" +
    				" for element-wise multiplication.");
    	}
    	
    	int numRows = matrixA.getRowDimension();
    	int numCols = matrixA.getColumnDimension();
    	RealMatrix product = new Array2DRowRealMatrix( numRows, numCols );
    	double aVal, bVal;
    	
    	//loop through each row
    	for(int r = 0; r < numRows; r++){
    	    //multiply each element of A by same element of B
    		for( int c = 0; c < numCols; c++){
    			aVal = matrixA.getEntry(r, c);
    			bVal = matrixB.getEntry(r, c);	
    			product.setEntry(r, c, aVal*bVal );
    	    }
    	}
    	return product;
    }
    

    /**
     * This method will return vech(matrix).
     * matrix must be symmetric in order to perform vech operation.
     * @param RealMatrix matrix
     * @return vech(matrix)
     */
	public static RealMatrix getVechMatrix(RealMatrix matrix)
	throws IllegalArgumentException{
    	if( matrix == null ||! isSymmetric( matrix ) ){
    		throw new IllegalArgumentException("Matrix must be non-null and " +
    				"symmetrical.");
    	}
    	int newRow = 0;
		int numRows = matrix.getRowDimension();
		RealMatrix vech = new Array2DRowRealMatrix(numRows*(numRows+1)/2, 1);
    	for( int c = 0; c < matrix.getColumnDimension();c++){
    		for(int r = c; r < matrix.getRowDimension(); r++, newRow++){
    			vech.setEntry(newRow, 0, matrix.getEntry(r, c));
    		}
    	}
    	return vech;
    }
    
    /**
     * Convenience method which will return a boolean which indicates whether
     * or not the dimensions of the two matrices are equal.  i.e., the method
     * will return true if the number of rows in matrixA is equal to the 
     * number of rows in matrixB, and the number of columns in matrixA equals
     * the number of columns in matrixB.
     * @param matrixA
     * @param matrixB
     * @return boolean indicating whether or not the matrix dimensions are equal.
     */
    public static boolean areDimensionsEqual(RealMatrix matrixA, RealMatrix matrixB){
    	int numRowsA = matrixA.getRowDimension();
    	int numColsA = matrixA.getColumnDimension();
    	int numRowsB = matrixB.getRowDimension();
    	int numColsB = matrixB.getColumnDimension();
    	return (numRowsA == numRowsB && numColsA == numColsB);
    }
    
	/**
	 * The method determines if the given matrix is symmetric, i.e.,
	 * if (r,c) == (c,r) in every case.
	 * @param RealMatrix matrix
	 * @return true if the matrix is symmetric, false if not symmetric.
	 */
    public static boolean isSymmetric(RealMatrix matrix){
    	int numRows = matrix.getRowDimension();
    	int numCols = matrix.getColumnDimension();
    	
    	//if not square, can't be symmetrical
    	if( numRows != numCols) return false;
    	
    	for( int r = 0; r < numRows; r++){
    		for( int c = 0; c < numCols; c++){
    			//test for symmetry
    			if( matrix.getEntry(r, c) != matrix.getEntry(c, r)){
    				return false;
    			}
    		}
    	}
    	return true;
    }
    
    /**
     * The method determines if the given matrix is positive definite.
     * The matrix must be square.
     * @param matrix
     * @param eigenTolerance is a double.  @see MatrixConstants.EIGEN_TOLERANCE
     * @return true if the matrix is positive definite.
     */
    public static boolean isPositiveDefinite(RealMatrix matrix, double eigenTolerance){
    	if( matrix == null || ! matrix.isSquare()){
    		throw new IllegalArgumentException("Matrix must be non-null, " +
    				"square. ");
    	}
    	double[] eigenValues = new EigenDecompositionImpl( matrix, eigenTolerance)
    	.getRealEigenvalues();

        // if all eigenValues are positive, we return true
        boolean isPositiveDefinite = true;
        for( int i = 0; i < eigenValues.length; i++){
        	if( eigenValues[i] <= 0){
        		isPositiveDefinite = false;
        		break;
        	}
        }
        return isPositiveDefinite;
    }
}
