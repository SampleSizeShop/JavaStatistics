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

import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;

/**
 * Class describing an essence design matrix, which is an abbreviated representation
 * of the full design matrix for a general linear multivariate model.
 * <p>
 * The EssenceMatrix class includes the following meta information used to
 * generate the full design matrix:
 * <ul>
 * <li>Column meta data: indicates if a predictor is fixed or random
 * <li>Row meta data: indicates the ratio of group sizes
 * </ul>
 * 
 * @author Sarah Kreidler
 *
 */
public class EssenceMatrix extends Array2DRowRealMatrix
{
    static final long serialVersionUID = 12357911L;
    
    // Indicates if the specified predictor is fixed or random.
    ColumnMetaData[] columnMetaData = null;
    // indicates how many times the row should be repeated in the full design matrix
    RowMetaData[] rowMetaData = null;
    
    // count of #random predictors
    int randomPredictorCount = 0;
    
    // random seed for expanding random covariates in the essence matrix
    int randomSeed = 1234;
    
    // per group sample size 
    int groupSampleSize = 1;
    
    /**
     * Constructor.  Creates an essence matrix from a real matrix.
     * Column and row meta data are set to default values
     * 
     * @param matrix RealMatrix containing essence matrix data
     */
    public EssenceMatrix(RealMatrix matrix)
    {
        this(matrix.getData());
    }
    
    /**
     * Constructor.  Creates an essence matrix with the specified data.
     * Column and row meta data are set to default values
     * 
     * @param data 2D array containing matrix data
     */
    public EssenceMatrix(double[][] data)
    {
        super(data);
        
        int columns = this.getColumnDimension();
        int rows = this.getRowDimension();
        
        // create default values for the meta information
        if (columns > 0)
        {
            columnMetaData = new ColumnMetaData[columns];
            for(int i = 0; i < columnMetaData.length; i++) 
                columnMetaData[i] = new ColumnMetaData();
        }
        if (rows > 0)
        {
            rowMetaData = new RowMetaData[rows];
            for(int i = 0; i < rowMetaData.length; i++) 
                rowMetaData[i] = new RowMetaData(1,1);
        }
    }
    
    /**
     * Constructor.  Creates an essence matrix with specified data and
     * meta information for predictors and rows.
     * <p>
     * The number of entries in the row and column meta data must match
     * the rows and columns in the data, respectively
     * 
     * @param data matrix contents
     * @param predictorTypes types of predictor variables in each column
     * @param rowRepetitions
     */
    public EssenceMatrix(double[][] data, ColumnMetaData[] columnMetaData,
            RowMetaData[] rowMetaData)
    throws IllegalArgumentException
    {
        super(data);
        
        if (columnMetaData != null)
        {
            this.columnMetaData = columnMetaData;
            for(ColumnMetaData colMDEntry: columnMetaData)
            {
                if (colMDEntry.getPredictorType() == PredictorType.RANDOM)
                    randomPredictorCount++;
            }
        }
        if (rowMetaData != null)
        {
            this.rowMetaData = rowMetaData;
        }
    }
        
    /**
     * Returns the number of random predictors in the model
     * 
     * @return number of random predictors in the model
     */
    public int getRandomPredictorCount()
    {
        return randomPredictorCount;
    }
    
    /**
     * Set the per group sample size.  If the group sizes are not equal
     * then the specified sample size will be multiplied by the group size
     * ratio specified in the row meta data for the design matrix.
     * <p/>
     * For example, for a 3x3 design matrix with group sizes 1:2:1 and
     * a groupN of 10, the actual group sizes will be 10, 20, and 10.     * 
     * 
     * @param groupN
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
     * @returns per group sample size
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
        int fullColumns = getColumnDimension();
        Array2DRowRealMatrix fullDesign = 
            new Array2DRowRealMatrix(fullRows, fullColumns);

        // copy in the columns
        // case 1: if predictor in a given column is fixed, the values are copied from 
        // the design matrix
        // case 2: if the predictor is random, the values are set to a random value from 
        // a normal curve with the mean/variance specified in the column meta data
        for(int col = 0; col < fullColumns; col++)
        {
            fillColumn(col, fullDesign, groupSampleSize);
        }
        return fullDesign;
    }
    
    /**
     * Expands the essence matrix into full design matrix for sample size
     * calculations.  
     * 
     * Assumes that any row meta data indicates the ratio of unique rows
     * in the full design (used to allow unequal group sizes in sample size
     * estimation)
     * 
     * @param totalN
     * @return
     * @throws IllegalArgumentException
     */
    public RealMatrix getFullDesignMatrix(int totalN)
    throws IllegalArgumentException
    {
        // allocate the full design matrix
        // #rows = total of repetitions for each unique row
        // #columns = number of columns in design matrix
        int repMultiplier = getRepetitionMultiplier(totalN);
        int fullRows = getMinimumSampleSize() * repMultiplier;
        int fullColumns = getColumnDimension();
        Array2DRowRealMatrix fullDesign = 
            new Array2DRowRealMatrix(fullRows, fullColumns);

        // copy in the columns
        // case 1: if predictor in a given column is fixed, the values are copied from 
        // the design matrix
        // case 2: if the predictor is random, the values are set to a random value from 
        // a normal curve with the mean/variance specified in the column meta data
        for(int col = 0; col < fullColumns; col++)
        {
            fillColumn(col, fullDesign, repMultiplier);
        }
        return fullDesign;
    }
    
    /**
     * Fills in a single column in the full design matrix
     * 
     * @param column
     * @param fullDesign
     * @param repMultiplier
     */
    private void fillColumn(int column, RealMatrix fullDesign, int repMultiplier)
    {
        ColumnMetaData colMD = columnMetaData[column];
        
        // if the column represents a random predictor, build a normal distribution
        // from which to pull random values
        Normal dist = null;
        if (colMD.getPredictorType() == PredictorType.RANDOM)
        {
            // note, the jsc library takes a standard deviation, not a variance so
            // we take the square root
            dist = new Normal(colMD.getMean(), Math.sqrt(colMD.getVariance()));
            dist.setSeed(randomSeed);
        }
        
        int essenceRow = 0;
        int reps = repMultiplier * rowMetaData[essenceRow].getRatio();
        for(int row = 0; row < fullDesign.getRowDimension(); row++)
        {
            // check if we need to move on to the next row in the essence matrix
            if (reps <= 0) 
            {
                essenceRow++;
                reps = repMultiplier * rowMetaData[essenceRow].getRatio();
            }
                
            // fill in the data
            if (dist == null)
            {
                // fixed predictor
                fullDesign.setEntry(row, column, this.getEntry(essenceRow, column));
            }
            else
            {
                // random predictor
                fullDesign.setEntry(row, column, dist.random());
            }
            
            // decrement the number of reps remain for this row
            reps--;    
        }
    }
    
    /**
     * Set the meta data for a specific column.  Throws an illegal argument
     * exception if the column index is out of bounds.
     * 
     * @param column column index (0 = first column)
     * @param metaData column meta data object
     * @throws IllegalArgumentException
     */
    public void setColumnMetaData(int column, ColumnMetaData metaData)
    throws IllegalArgumentException
    {
        if (column < 0 || column >= columnMetaData.length)
            throw new IllegalArgumentException("Requested column [" + column + "] is outside matrix bounds");
        else
        {
            PredictorType currentType = columnMetaData[column].getPredictorType();
            // adjust the random predictor count
            if (currentType == PredictorType.FIXED && 
                    metaData.getPredictorType() == PredictorType.RANDOM)
                randomPredictorCount++;
            else if (currentType == PredictorType.RANDOM && 
                    metaData.getPredictorType() == PredictorType.FIXED)
                randomPredictorCount--;
            // reset the meta data
            columnMetaData[column] = metaData;
        }
    }
    
    /**
     * Set the meta data for all columns.  
     * 
     * @param column column index (0 = first column)
     * @param metaData column meta data object
     * @throws IllegalArgumentException
     */
    public void setColumnMetaData(ColumnMetaData[] metaData)
    throws IllegalArgumentException
    {
        if (metaData == null || metaData.length != this.getColumnDimension())
            throw new IllegalArgumentException("Invalid column meta data.  Must have same number of columns as matrix data");
        columnMetaData = metaData;
        // reset the random predictor count
        randomPredictorCount = 0;
        for(ColumnMetaData colMDEntry: columnMetaData)
        {
            if (colMDEntry.getPredictorType() == PredictorType.RANDOM)
                randomPredictorCount++;
        }
    }
    
    /**
     * Get the meta data for a specific column.  Throws an illegal argument
     * exception if the column index is out of bounds.
     * 
     * @param column column index (0 = first column)
     * @return meta data object for the column
     * @throws IllegalArgumentException
     */
    public ColumnMetaData getColumnMetaData(int column)
    throws IllegalArgumentException
    {
        if (column < 0 || column >= columnMetaData.length)
            throw new IllegalArgumentException("Requested column [" + column + "] is outside matrix bounds");
        else
            return columnMetaData[column];
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
     * @param row row index (0 = first row)
     * @param metaData row meta data object
     * @throws IllegalArgumentException
     */
    public void setRowMetaData(RowMetaData[] metaData)
    throws IllegalArgumentException
    {
        if (metaData == null || metaData.length != this.getRowDimension())
            throw new IllegalArgumentException("Invalid row meta data.  Must have same number of rows as matrix data");
        rowMetaData = metaData;
    }
    
    /**
     * Get the meta data for a specific row.  Throws an illegal argument
     * exception if the row index is out of bounds.
     * 
     * @param row row index (0 = first row)
     * @return meta data object for the row
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
     * Determine the multiplier for the row repetitions when the essence
     * matrix specifies an unequal ratio of group sizes.
     *  
     * @param totalN total rows in full design matrix
     * @return repetition multiplier for each unique row in essence matrix
     */
    private int getRepetitionMultiplier(int totalN)
    {
        int minSampleSize = getMinimumSampleSize();
        if (totalN < minSampleSize)
        {
            return 1;
        }
        else
        {
            int remainder = totalN % minSampleSize;          
            if (remainder != 0) totalN += minSampleSize - remainder;
            return (int) Math.round((double) totalN / (double) minSampleSize);

        }
    }
    
	public int getRandomSeed()
	{
		return randomSeed;
	}

	public void setRandomSeed(int randomSeed)
	{
		this.randomSeed = randomSeed;
	}
    
	public RealMatrix getFullDesignFixed()
	{
	    int designRows = getFullDesignMatrix().getRowDimension();
	    int[] rows = new int[designRows];
	    int[] cols = new int[getColumnDimension() - getRandomPredictorCount()];
	    // include all rows
	    for(int r = 0; r < designRows; r++) rows[r] = r;
	    // only include columns for fixed predictors
	    int i = 0;
	    for(int c = 0; c < getColumnDimension(); c++)
	    {
	        if (columnMetaData[c].getPredictorType() == PredictorType.FIXED)
	        {
	            cols[i] = c;
	            i++;
	        }
	    }
	    return getFullDesignMatrix().getSubMatrix(rows, cols); 
	}

}
