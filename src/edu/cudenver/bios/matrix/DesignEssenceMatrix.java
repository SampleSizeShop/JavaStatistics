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
public class DesignEssenceMatrix extends FixedRandomMatrix
{
    static final long serialVersionUID = 12357911L;
    
    RealMatrix fixedMatrix = null;
    RealMatrix randomMatrix = null;
    
    // indicates how many times the row should be repeated in the full design matrix
    RowMetaData[] rowMetaData = null;
    
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
    public DesignEssenceMatrix(double[][] fixedData, double[][] randomData,
            RowMetaData[] rowMetaData)
    throws IllegalArgumentException
    {
        super(fixedData, randomData, true);
        if (randomMatrix != null && randomMatrix.getColumnDimension() > 1)
            throw new IllegalArgumentException("Only a single random covariate is allowed");
        if (combinedMatrix.getRowDimension() != rowMetaData.length)
            throw new IllegalArgumentException("Row meta data array does not match the number of rows in the data");
        this.rowMetaData = rowMetaData;
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
        int fullColumns = combinedMatrix.getColumnDimension();
        Array2DRowRealMatrix fullDesign = 
            new Array2DRowRealMatrix(fullRows, fullColumns);

        // copy in the columns
        // case 1: if predictor in a given column is fixed, the values are copied from 
        // the design matrix
        // case 2: if the predictor is random, the values are set to a random value from 
        // a normal curve with the mean/variance specified in the column meta data
        for(int col = 0; col < fullColumns; col++)
        {
            fillColumn(col, fullDesign);
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
    private void fillColumn(int column, RealMatrix fullDesign)
    {
//
//        
//        // if the column represents a random predictor, build a normal distribution
//        // from which to pull random values
//        Normal dist = null;
//        if (colMD.getPredictorType() == PredictorType.RANDOM)
//        {
//            // note, the jsc library takes a standard deviation, not a variance so
//            // we take the square root
//            dist = new Normal(colMD.getMean(), Math.sqrt(colMD.getVariance()));
//            dist.setSeed(randomSeed);
//        }
//        
//        int essenceRow = 0;
//        int reps = repMultiplier * rowMetaData[essenceRow].getRatio();
//        for(int row = 0; row < fullDesign.getRowDimension(); row++)
//        {
//            // check if we need to move on to the next row in the essence matrix
//            if (reps <= 0) 
//            {
//                essenceRow++;
//                reps = repMultiplier * rowMetaData[essenceRow].getRatio();
//            }
//                
//            // fill in the data
//            if (dist == null)
//            {
//                // fixed predictor
//                fullDesign.setEntry(row, column, this.getEntry(essenceRow, column));
//            }
//            else
//            {
//                // random predictor
//                fullDesign.setEntry(row, column, dist.random());
//            }
//            
//            // decrement the number of reps remain for this row
//            reps--;    
//        }
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
        if (metaData == null || metaData.length != combinedMatrix.getRowDimension())
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
//    
//    /**
//     * Determine the multiplier for the row repetitions when the essence
//     * matrix specifies an unequal ratio of group sizes.
//     *  
//     * @param totalN total rows in full design matrix
//     * @return repetition multiplier for each unique row in essence matrix
//     */
//    private int getRepetitionMultiplier(int totalN)
//    {
//        int minSampleSize = getMinimumSampleSize();
//        if (totalN < minSampleSize)
//        {
//            return 1;
//        }
//        else
//        {
//            int remainder = totalN % minSampleSize;          
//            if (remainder != 0) totalN += minSampleSize - remainder;
//            return (int) Math.round((double) totalN / (double) minSampleSize);
//
//        }
//    }
    
    public int getRandomSeed()
    {
        return randomSeed;
    }

    public void setRandomSeed(int randomSeed)
    {
        this.randomSeed = randomSeed;
    }
    


}
