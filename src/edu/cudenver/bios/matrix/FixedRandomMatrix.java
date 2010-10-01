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

/**
 * Matrix which contains fixed and random components.  The combined
 * matrix may be produced by concatenating the fixed and random
 * submatrices either vertically or horizontally.
 * 
 * @author Sarah Kreidler
 *
 */
public class FixedRandomMatrix
{
    protected RealMatrix fixedMatrix;
    protected RealMatrix randomMatrix;
    protected RealMatrix combinedMatrix;
    protected boolean combineHorizontal = true;
    
    /**
     * Constructor, create a fixed/random matrix from the specified data arrays
     * 
     * @param fixedData submatrix of fixed predictor information
     * @param randomData submatrix of random predictor information
     * @param combineHorizontal indicates if the combined matrix should be
     * formed by horizontal or vertical concatenation
     * @throws IllegalArgumentException
     */
    public FixedRandomMatrix(double[][] fixedData, double[][] randomData,
            boolean combineHorizontal) throws IllegalArgumentException
    {
        if (fixedData == null && randomData == null)
            throw new IllegalArgumentException("No data specified");
        
        this.combineHorizontal = combineHorizontal;
   
        // determine the number of rows/columns in the combination of the two matrixes
        int totalRows = 0;
        int totalColumns = 0;
        if (fixedData != null) 
        {
        	this.fixedMatrix = new Array2DRowRealMatrix(fixedData);
        	totalRows += fixedMatrix.getRowDimension();
        	totalColumns = fixedMatrix.getColumnDimension();
        }
        if (randomData != null)
        {
        	this.randomMatrix = new Array2DRowRealMatrix(randomData);
        	if (combineHorizontal)
        	{
        		totalColumns += randomMatrix.getColumnDimension();
        		if (fixedData == null) totalRows += randomMatrix.getRowDimension();
        	}
        	else
        	{
        		totalRows += randomMatrix.getRowDimension();
        		if (fixedData == null) totalColumns += randomMatrix.getColumnDimension();
        	}
        }

        // before we continue, make sure the matrices can be combined as specified
        if (fixedMatrix != null && randomMatrix != null)
        {
        	if (combineHorizontal && fixedMatrix.getRowDimension() != randomMatrix.getRowDimension())
        	{
                throw new IllegalArgumentException("Fixed and random matrices must have same number of rows for horizontal combination");
        	}
        	else if (!combineHorizontal && fixedMatrix.getColumnDimension() != randomMatrix.getColumnDimension())
        	{
                throw new IllegalArgumentException("Fixed and random matrices must have same number of columns for vertical combination");
        	}
        }
        
        // allocate the combined matrix
        combinedMatrix = new Array2DRowRealMatrix(totalRows, totalColumns);
        // fill in the data
        if (fixedMatrix != null)
        {
        	combinedMatrix.setSubMatrix(fixedData, 0, 0);
        }
        if (randomMatrix != null)
        {
        	if (fixedMatrix == null)
        	{
        		combinedMatrix.setSubMatrix(randomData, 0, 0);
        	}
        	else
        	{
        		if (combineHorizontal)
        			combinedMatrix.setSubMatrix(randomData, 0, fixedMatrix.getColumnDimension());
        		else
        			combinedMatrix.setSubMatrix(randomData, fixedMatrix.getRowDimension(), 0);        			
        	}
        }
    }
    
    /**
     * Get the fixed submatrix
     * @return fixed submatrix
     */
    public RealMatrix getFixedMatrix()
    {
        return fixedMatrix;
    }
    
    /**
     * Get the random submatrix
     * @return random submatrix
     */
    public RealMatrix getRandomMatrix()
    {
        return randomMatrix;
    }

    /**
     * Get the matrix resulting from the concatenation
     * of the fixed and random submatrices
     * @return combined matrix
     */
    public RealMatrix getCombinedMatrix()
    {
        return combinedMatrix;
    }
    
    /**
     * Determine if this fixed/random matrix includes a fixed component
     * @return true if fixed component is present
     */
    public boolean hasFixed()
    {
        return (fixedMatrix != null);
    }
    
    /**
     * Determine if this fixed/random matrix includes a random component
     * @return true if random component is present
     */
    public boolean hasRandom()
    {
        return (randomMatrix != null);
    }
    
    /**
     * Perform scalar multiplication on a fixed / random matrix and
     * return the resulting scaled matrix.  The scale factor can be applied 
     * to the entire matrix, or just the fixed submatrix.
     * 
     * @param scale scale factor
     * @param fixedOnly if true, only the fixed submatrix will be scaled
     * @return scaled combined matrix
     */
    public RealMatrix scalarMultiply(double scale, boolean fixedOnly)
    {
    	if (fixedOnly)
    	{
    		RealMatrix fixedScaled = fixedMatrix.scalarMultiply(scale);
    		// update the combined matrix and return it
    		combinedMatrix.setSubMatrix(fixedScaled.getData(), 0, 0);
    		return combinedMatrix;
    	}
    	else
    	{
    		return combinedMatrix.scalarMultiply(scale);
    	}
    }
    
    /**
     * Update the values in the random submatrix.  The matrix of new values
     * must match the dimensions of the random submatrix.
     * 
     * @param updatedMatrix new values for the random submatrix
     */
    public void updateRandomMatrix(RealMatrix updatedMatrix)
    {
    	if (randomMatrix != null && 
    			updatedMatrix.getColumnDimension() == randomMatrix.getColumnDimension() &&
    			updatedMatrix.getRowDimension() == randomMatrix.getRowDimension())
    	{
    		double[][] data = updatedMatrix.getData();
    		randomMatrix.setSubMatrix(data, 0, 0);
    		int startRow = 0;
    		int startCol = 0;
    		if (fixedMatrix != null)
    		{
    			if (combineHorizontal)
    			{
    				startRow = 0;
    				startCol = fixedMatrix.getColumnDimension();
    			}
    			else
    			{
    				startRow = fixedMatrix.getRowDimension();
    				startCol = 0;
    			}
    		}
    		combinedMatrix.setSubMatrix(data, startRow, startCol);
    	}
    }
}
