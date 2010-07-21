package edu.cudenver.bios.matrix;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

public class FixedRandomMatrix
{
    protected RealMatrix fixedMatrix;
    protected RealMatrix randomMatrix;
    protected RealMatrix combinedMatrix;
    protected boolean combineHorizontal = true;
    
    public FixedRandomMatrix(double[][] fixedData, double[][] randomData,
            boolean combineHorizontal) throws IllegalArgumentException
    {
        if (fixedData == null && randomData == null)
            throw new IllegalArgumentException("No data specified");
        
        if (fixedData != null) this.fixedMatrix = new Array2DRowRealMatrix(fixedData);
        if (randomData != null) this.randomMatrix = new Array2DRowRealMatrix(randomData);
        this.combineHorizontal = combineHorizontal;
        
        double[][] combinedData = null;
        if (fixedData != null && randomData != null)
        {
            if (combineHorizontal)
            {
                if (this.fixedMatrix.getRowDimension() != this.randomMatrix.getRowDimension())
                    throw new IllegalArgumentException("fixed and random data must have the same number of rows");
                int rows = fixedMatrix.getRowDimension();
                int columns = fixedMatrix.getColumnDimension()+randomMatrix.getColumnDimension();
                combinedData = 
                    new double[rows][columns];
                for(int r = 0; r < rows; r++)
                {
                    int c = 0;
                    for(; c < fixedMatrix.getColumnDimension(); c++)
                    {
                        combinedData[r][c] = fixedData[r][c];
                    }
                    for(int randC = 0; randC < randomMatrix.getColumnDimension(); randC++, c++)
                    {
                        combinedData[r][c] = randomData[r][randC];
                    }
                }
                
            }
            else
            {
                if (this.fixedMatrix.getColumnDimension() != this.randomMatrix.getColumnDimension())
                    throw new IllegalArgumentException("fixed and random data must have the same number of columns");
                int rows = fixedMatrix.getRowDimension()+randomMatrix.getRowDimension();
                int columns = fixedMatrix.getColumnDimension();
                combinedData = 
                    new double[rows][columns];
                for(int c = 0; c < columns; c++)
                {
                    int r = 0;
                    for(; r < fixedMatrix.getRowDimension(); r++)
                    {
                        combinedData[r][c] = fixedData[r][c];
                    }
                    for(int randR = 0; randR < randomMatrix.getRowDimension(); randR++, r++)
                    {
                        combinedData[r][c] = randomData[r][randR];
                    }
                }
            }
        }
        else if (fixedData != null)
        {
            combinedData = fixedData;
        }
        else
        {
            combinedData = randomData;
        }
        
        combinedMatrix = new Array2DRowRealMatrix(combinedData);
    }
    
    public RealMatrix getFixedMatrix()
    {
        return fixedMatrix;
    }
    
    public RealMatrix getRandomMatrix()
    {
        return randomMatrix;
    }

    public RealMatrix getCombinedMatrix()
    {
        return combinedMatrix;
    }
}
