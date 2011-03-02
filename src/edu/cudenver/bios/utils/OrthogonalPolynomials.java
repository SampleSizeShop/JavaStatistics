package edu.cudenver.bios.utils;

import java.text.DecimalFormat;
import java.util.HashSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.QRDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;

public class OrthogonalPolynomials
{

	/**
	 * Computes orthogonal polynomial contrasts for the specified data values.  Currently only
	 * supports fitting (not prediction contrasts).  Based on the poly command from R 
	 * (http://cran.r-project.org/)
	 * 
	 * 	@param x the points at which the polynomials will be evaluated
	 * @param maxDegree contrasts will be computed for degrees 1 to maxDegree
	 * @return matrix containing 1st, 2nd,...maxDegree-th degree contrasts in each column
	 * @throws IllegalArgumentException
	 */
	public static RealMatrix orthogonalPolynomialCoefficients(double[] x, int maxDegree)
	throws IllegalArgumentException
	{
		if (x == null) throw new IllegalArgumentException("no data specified");
		if (maxDegree < 1)
			throw new IllegalArgumentException("max polynomial degree must be greater than 1");
		// count number of unique values
		HashSet<Double> s = new HashSet<Double>();
		for (double i : x) s.add(i);
		int uniqueCount = s.size();
		if (maxDegree >= uniqueCount)
			throw new IllegalArgumentException("max polynomial degree must be less than the number of unique points");
		
		// center the data
		double xBar = StatUtils.mean(x);
		double[] xCentered = new double[x.length];
		for(int i = 0; i < x.length; i++) xCentered[i] = x[i] - xBar;
		// compute an "outer product" of the centered x vector and a vector 
		// containing the sequence 0 to maxDegree-1, but raise the x values
		// to the power in the sequence array
		double[][] xOuter = new double[x.length][maxDegree+1];
		int row = 0;
		for(double xValue: xCentered)
		{
			for(int col = 0; col <= maxDegree; col++)
			{
				xOuter[row][col] = Math.pow(xValue, col);
			}
			row++;
		}		
		// do some mysterious QR decomposition stuff.  See Emerson (1968)
		RealMatrix outerVector = new Array2DRowRealMatrix(xOuter);
		QRDecompositionImpl qrDecomp = new QRDecompositionImpl(outerVector);
		double[][] zData = new double[qrDecomp.getR().getRowDimension()][qrDecomp.getR().getColumnDimension()]; 
		for(row = 0; row < qrDecomp.getR().getRowDimension(); row++)
		{
			for(int col = 0; col < qrDecomp.getR().getColumnDimension(); col++)
			{
				if (row == col) 
					zData[row][col] = qrDecomp.getR().getEntry(row, col);
				else
					zData[row][col] = 0;
			}
		}
		RealMatrix z = new Array2DRowRealMatrix(zData);
		RealMatrix raw = qrDecomp.getQ().multiply(z);
		
		// column sum of squared elements in raw
		double[] normalizingConstants = new double[raw.getColumnDimension()];
		for(int col = 0; col < raw.getColumnDimension(); col++)
		{
			normalizingConstants[col] = 0;
			for(row = 0; row < raw.getRowDimension(); row++)
			{
				double value = raw.getEntry(row, col);
				normalizingConstants[col] += value*value;
			}
		}
		
		// now normalize the raw values
		for(int col = 0; col < raw.getColumnDimension(); col++)
		{
			double normalConstantSqrt = Math.sqrt(normalizingConstants[col]);
			for(row = 0; row < raw.getRowDimension(); row++)
			{
				raw.setEntry(row, col, raw.getEntry(row, col) / normalConstantSqrt);
			}
		}

		return raw;
	}
	
	public static RealMatrix withinSubjectContrastOneFactor()
	{
		Array2DRowRealMatrix orpolMatrix = null;
		
		return orpolMatrix;
	}
	
	public static RealMatrix withinSubjectContrastTwoFactor()
	{
		Array2DRowRealMatrix orpolMatrix = null;
		
		return orpolMatrix;
	}
	
	public static RealMatrix withinSubjectContrastThreeFactor()
	{
		Array2DRowRealMatrix orpolMatrix = null;
		
		return orpolMatrix;
	}
	
	/**
	 * Mini test program for the orthogonal polynomial coefficients
	 * TODO: move this to unit test
	 * @param args null
	 */
	public static void main(String[] args)
	{
		double[] x = {1,2,3,4,5,6};
		RealMatrix poly = OrthogonalPolynomials.orthogonalPolynomialCoefficients(x,4);
	
	    DecimalFormat Number = new DecimalFormat("#0.000");
		for(int row = 0; row < poly.getRowDimension(); row++)
		{
			for(int col= 0; col < poly.getColumnDimension(); col++)
			{
				System.out.print(Number.format(poly.getEntry(row, col)) + "\t");
			}
			System.out.print("\n");
		}
	}
	
}
