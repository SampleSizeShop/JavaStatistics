package edu.cudenver.bios.matrix;

import jsc.distributions.Normal;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

/**
 * Matrix of random normals with the specified covariance matrix
 */
public class RandomErrorMatrix
{
	private static final double RELATIVE_SYMMETRY_THRESHOLD = 
	    CholeskyDecompositionImpl.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;
	private static final double POSITIVITY_THRESHOLD = 
	    CholeskyDecompositionImpl.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
	
	protected long seed = 1234;
	protected RealMatrix matrix;
	protected RealMatrix sigma;
	protected Normal normalDist;

	public RandomErrorMatrix(int rows, int cols, RealMatrix sigma)
	{
		matrix = new Array2DRowRealMatrix(rows, cols);
		this.sigma = sigma;
	}

	public void setSeed(long seed)
	{
		this.seed = seed;
		normalDist.setSeed(seed);
	}
	
    /**
     * Simulate the error matrix in the Y = X * beta + e
     * 
     * @param normalDist normal distribution object for generating random samples
     * @param error matrix to hold random normal samples
     * @param rows number of rows in the random normal samples matrix
     * @param columns number of columns in the random normal samples matrix
     * @param sigma the covariance matrix for the error term
     * @return a random instance of the 'e' matrix in the model
     */
	public RealMatrix random()
	throws IllegalArgumentException
	{
        // build a matrix of random values from a standard normal
        // the number of rows = #subjects (rows) in the full design matrix
        // the number of columns = #outcome variables (i.e. columns in beta)
        for(int rowIndex = 0; rowIndex < matrix.getRowDimension(); rowIndex++)
        {
            for(int columnIndex = 0; columnIndex < matrix.getColumnDimension(); columnIndex++)
            {
            	matrix.setEntry(rowIndex, columnIndex, normalDist.random()); 
            }
        }
        
        // take the square root of the sigma matrix via cholesky decomposition
        try
        {            
            RealMatrix sqrtMatrix = 
                new CholeskyDecompositionImpl(sigma, 
                        RELATIVE_SYMMETRY_THRESHOLD,
                        POSITIVITY_THRESHOLD).getLT();
            return matrix.multiply(sqrtMatrix); 
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
	}

}
