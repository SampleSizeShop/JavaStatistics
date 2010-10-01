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

import jsc.distributions.Normal;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

/**
 * Class describing an essence design matrix, which is an abbreviated representation
 * of the full design matrix for a general linear multivariate model.
 * <p>
 * The EssenceMatrix class includes the following meta information used to
 * generate the full design matrix:
 * <ul>
 * <li>Random column meta data: for columns representing predictors with normally distributed
 * values, this meta data indicates the mean and variance of the distribution
 * <li>Row meta data: indicates the ratio of group sizes for each unique row in the design matrix
 * </ul>
 * 
 * @author Sarah Kreidler
 *
 */
public class DesignEssenceMatrix extends FixedRandomMatrix
{
    static final long serialVersionUID = 12357911L;
    
    // indicates how many times the row should be repeated in the full design matrix
    RowMetaData[] rowMetaData = null;
    
    // contains means/variance for each random predictor
    RandomColumnMetaData[] randomColMetaData = null;
    
    // random seed for expanding random covariates in the essence matrix
    int randomSeed = 1234;
    
    // per group sample size 
    int groupSampleSize = 1;
    
    /**
     * Constructor.  Creates an essence matrix from data sets for
     * submatrices representing fixed predictors and random predictors.
     * If random predictors are present, then distribution information should
     * be specified in the randomColMetaData array.
     * 
     * @param fixedData portion of design matrix representing fixed predictors
     * @param rowMetaData relative group sizes for each unique row in the design matrix
     * @param randomData portion of design matrix representing random covariates
     * @param randomColMetaData distribution (mean, variance) information for the random covariates
     */
    public DesignEssenceMatrix(double[][] fixedData, RowMetaData[] rowMetaData, 
    		double[][] randomData, RandomColumnMetaData[] randomColMetaData)
    throws IllegalArgumentException
    {
        super(fixedData, randomData, true);
        if (randomMatrix != null && randomMatrix.getColumnDimension() > 1)
            throw new IllegalArgumentException("Only a single random covariate is allowed");
        if (combinedMatrix.getRowDimension() != rowMetaData.length)
            throw new IllegalArgumentException("Row meta data array does not match the number of rows in the data");
        this.rowMetaData = rowMetaData;
        this.randomColMetaData = randomColMetaData;
    }

    /**
     * Set the per group sample size.  If the group sizes are not equal
     * then the specified sample size will be multiplied by the group size
     * ratio specified in the row meta data for the design matrix.
     * <p/>
     * For example, for a 3x3 design matrix with group sizes 1:2:1 and
     * a groupN of 10, the actual group sizes will be 10, 20, and 10.     * 
     * 
     * @param groupSampleSize new per group sample size
     */
    public void setGroupSampleSize(int groupSampleSize)
    throws IllegalArgumentException
    {
        if (groupSampleSize <= 0) 
            throw new IllegalArgumentException("Per group sample size must be positive");
        
        this.groupSampleSize = groupSampleSize;
    }
    
    /**
     * Get the per group sample size.  For non-equal group sizes, this
     * function will return the size of the smallest group.
     * 
     * @return per group sample size
     */
    public int getGroupSampleSize()
    {
        return groupSampleSize;
    }
    
    /**
     * Expands the essence matrix into full design matrix for power 
     * calculations.  
     * 
     * Assumes that row meta data indicates the actual number
     * of repetitions of that row in the full design
     * 
     * @return design matrix
     * @throws IllegalArgumentException
     */
    public RealMatrix getFullDesignMatrix()
    {
        // allocate the full design matrix
        // #rows = total of repetitions for each unique row
        // #columns = number of columns in design matrix
        int fullRows = getTotalSampleSize();
        int fullColumns = combinedMatrix.getColumnDimension();
        Array2DRowRealMatrix fullDesignMatrix = 
            new Array2DRowRealMatrix(fullRows, fullColumns);

        // copy in the columns
        // case 1: if predictor in a given column is fixed, the values are copied from 
        // the design matrix
        // case 2: if the predictor is random, the values are set to a random value from 
        // a normal curve with the mean/variance specified in the column meta data
        int fullDesignColumn = 0;
        if (fixedMatrix != null)
        {
        	for(int col = 0; col < fixedMatrix.getColumnDimension(); col++, fullDesignColumn++)
        	{
        		fillFixedColumn(col, fullDesignColumn, fullDesignMatrix);
        	}
        }
        if (randomMatrix != null)
        {
        	for(int col = 0; col < randomMatrix.getColumnDimension(); col++, fullDesignColumn++)
        	{
        		fillRandomColumn(col, fullDesignColumn, fullDesignMatrix);
        	}
        }
        return fullDesignMatrix;
    }

    /**
     * Expands the essence matrix into full design matrix for power 
     * calculations using the fixed component only
     * 
     * Assumes that row meta data indicates the actual number
     * of repetitions of that row in the full design
     * 
     * @return design matrix
     * @throws IllegalArgumentException
     */
    public RealMatrix getFullDesignMatrixFixed()
    {
    	if (fixedMatrix == null) return null;
        // allocate the full design matrix
        // #rows = total of repetitions for each unique row
        // #columns = number of columns in design matrix
        int fullRows = getTotalSampleSize();
        int fullColumns = fixedMatrix.getColumnDimension();
        Array2DRowRealMatrix fullDesignMatrix = 
            new Array2DRowRealMatrix(fullRows, fullColumns);

        // copy in the columns
        // case 1: if predictor in a given column is fixed, the values are copied from 
        // the design matrix
        // case 2: if the predictor is random, the values are set to a random value from 
        // a normal curve with the mean/variance specified in the column meta data
        int fullDesignColumn = 0;
        for(int col = 0; col < fixedMatrix.getColumnDimension(); col++, fullDesignColumn++)
        {
            fillFixedColumn(col, fullDesignColumn, fullDesignMatrix);
        }
        return fullDesignMatrix;
    }
    
    /**
     * Update the random columns in the specified design matrix with a new realization
     * of random values
     * 
     * @param fullDesignMatrix an instance of the full design matrix
     */
    public void updateRandomColumns(RealMatrix fullDesignMatrix)
    throws IllegalArgumentException
    {
    	if (fullDesignMatrix.getColumnDimension() != combinedMatrix.getColumnDimension())
    		throw new IllegalArgumentException("Design matrix instance has invalid #columns");
    	int fullDesignColumn = 0;
    	if (fixedMatrix != null) fullDesignColumn = fixedMatrix.getColumnDimension();
    	
        if (randomMatrix != null)
        {
        	for(int col = 0; col < randomMatrix.getColumnDimension(); col++, fullDesignColumn++)
        	{
        		fillRandomColumn(col, fullDesignColumn, fullDesignMatrix);
        	}
        }
    	
    }
    
    /**
     * Fill a fixed column in the design matrix
     * 
     * @param fixedColumn column index in fixed submatrix
     * @param fullColumn column index in full design matrix
     * @param fullDesign full design matrix
     */
    private void fillFixedColumn(int fixedColumn, int fullColumn, RealMatrix fullDesign)
    {       
    	int essenceRow = 0;
    	int reps = groupSampleSize * rowMetaData[essenceRow].getRatio();
    	for(int row = 0; row < fullDesign.getRowDimension(); row++)
    	{
    		// check if we need to move on to the next row in the essence matrix
    		if (reps <= 0) 
    		{
    			essenceRow++;
    			reps = groupSampleSize * rowMetaData[essenceRow].getRatio();
    		}

    		// fill in the data
    		fullDesign.setEntry(row, fullColumn, fixedMatrix.getEntry(essenceRow, fixedColumn));
    		// decrement the number of reps remain for this row
    		reps--;    
    	}
    }
    
    /**
     * Fills in a random column in the full design matrix
     * 
     * @param randomColumn column index in random submatrix
     * @param fullColumn column index in full design matrix
     * @param fullDesign full design matrix
     */
    private void fillRandomColumn(int randomColumn, int fullColumn, RealMatrix fullDesign)
    {
        // if the column represents a random predictor, build a normal distribution
        // from which to pull random values
        Normal dist = null;
        // note, the jsc library takes a standard deviation, not a variance so
        // we take the square root
        dist = new Normal(randomColMetaData[randomColumn].getMean(), 
        		Math.sqrt(randomColMetaData[randomColumn].getVariance()));
        dist.setSeed(randomSeed);

        for(int row = 0; row < fullDesign.getRowDimension(); row++)
        {                
            // fill in the data
                fullDesign.setEntry(row, fullColumn, dist.random());
        }
    }
    
    /**
     * Set the meta data for a specific row.  Throws an illegal argument
     * exception if the row index is out of bounds.
     * 
     * @param row row index (0 = first row)
     * @param metaData row meta data object
     * @throws IllegalArgumentException
     */
    public void setRowMetaData(int row, RowMetaData metaData)
    throws IllegalArgumentException
    {
        if (row < 0 || row >= rowMetaData.length)
            throw new IllegalArgumentException("Requested row [" + row + "] is outside matrix bounds");
        else
            rowMetaData[row] = metaData;
    }
    
    /**
     * Set the meta data for all of the rows.  
     * 
     * @param metaData row meta data object
     * @throws IllegalArgumentException
     */
    public void setRowMetaData(RowMetaData[] metaData)
    throws IllegalArgumentException
    {
        if (metaData == null || metaData.length != combinedMatrix.getRowDimension())
            throw new IllegalArgumentException("Invalid row meta data.  Must have same number of rows as matrix data");
        rowMetaData = metaData;
    }
    
    /**
     * Get the meta data for a specific row.  Throws an illegal argument
     * exception if the row index is out of bounds.
     * 
     * @param row row index (0 = first row)
     * @throws IllegalArgumentException
     */
    public RowMetaData getRowMetaData(int row)
    throws IllegalArgumentException
    {
        if (row < 0 || row >= rowMetaData.length)
            throw new IllegalArgumentException("Requested row [" + row + "] is outside matrix bounds");
        else
            return rowMetaData[row];
    }    
    
    /**
     * Get the total number of rows in the full design matrix.
     * Calculated as the total number of repetitions specified
     * in the essence matrix.
     * 
     * @return total rows in full design matrix
     */
    public int getTotalSampleSize()
    {
        int count = 0;
        for(RowMetaData md : rowMetaData)
        {
            count += groupSampleSize * md.getRatio();
        }
        return count;
    }
    
    /**
     * Returns the minimum number of subjects
     * @return minimum sample size
     */
    public int getMinimumSampleSize()
    {
        int ratioCount = 0;
        for(RowMetaData md : rowMetaData)
        {
            ratioCount += md.getRatio();
        }
        return ratioCount;
    }
    
    /**
     * Get the seed for random generation of the Gaussian predictor values
     * @return seed
     */
    public int getRandomSeed()
    {
        return randomSeed;
    }

    /**
     * Set the seed for random generation of the Gaussian predictor values
     * @param randomSeed
     */
    public void setRandomSeed(int randomSeed)
    {
        this.randomSeed = randomSeed;
    }
    
    /**
     * Return the column meta data for the specified column in the random
     * matrix.
     * @param randomColumn column in the random matrix
     */
    public RandomColumnMetaData getColumnMetaData(int randomColumn)
    {
    	if (randomMatrix != null)
    	{
    		if (randomColumn >= 0 && randomColumn < randomMatrix.getColumnDimension())
    		{
    			return randomColMetaData[randomColumn];
    		}
    	}
    	return null;
    }

}
