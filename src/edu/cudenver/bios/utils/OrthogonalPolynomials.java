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
package edu.cudenver.bios.utils;

import java.util.HashSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.QRDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;

import edu.cudenver.bios.matrix.MatrixUtils;

/**
 * Class to generate orthogonal polynomial contrasts for unequally spaced
 * measurements.
 * 
 * @author Sarah Kreidler
 *
 */
public class OrthogonalPolynomials
{

	/**
	 * Computes orthogonal polynomial contrasts for the specified data values.  Currently only
	 * supports fitting (not prediction contrasts).  
	 * 
	 * 	@param x the points at which the polynomials will be evaluated
	 * @param maxDegree contrasts will be computed for degrees 1 to maxDegree
	 * @return matrix containing 0th,1st, 2nd,...maxDegree-th degree contrasts in each column
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

		RealMatrix z = MatrixUtils.toDiagonalMatrix(qrDecomp.getR());
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
	
	/**
	 * Create a within subject contrast (U) for polynomial trends
	 * in a single factor.
	 * 
	 * @param values values of factor
	 * @return polynomial contrast matrix
	 * @throws IllegalArgumentException
	 */
	public static RealMatrix withinSubjectContrast(double[] values)
	throws IllegalArgumentException
	{		
		if (values == null || values.length < 2)
			throw new IllegalArgumentException("must specify at least 2 values");

		int p = values.length;

		double mean = StatUtils.mean(values);
		double[] centered = new double[p];
		double scale = 0;
		for(int i = 0; i < p; i++) 
		{
			double v = values[i] - mean;
			centered[i] = v;
			scale += v*v;
		}
		for(int i = 0; i < p; i++) centered[i] /= Math.sqrt(scale);
		
		RealMatrix poly = OrthogonalPolynomials.orthogonalPolynomialCoefficients(centered,p-1);
		int[] rows = new int[poly.getRowDimension()];
		for(int i = 0; i < rows.length; i++) rows[i] = i;
		int[] cols = new int[poly.getColumnDimension()-1];
		for(int i = 0; i < cols.length; i++) cols[i] = i+1;
		return poly.getSubMatrix(rows, cols);	
	}
	
	/**
	 * Create a within subject contrast (U) for polynomial trends
	 * in two factors using a Kronecker product. It assumes Factor 1
	 *  varies slowly and that Factor 2 varies rapidly.    
	 * 
	 * @param factor1Values values of first factor
	 * @param factor2Values values of second factor
	 * @return 2 factor polynomial contrast matrix
	 * @throws IllegalArgumentException
	 */
	public static RealMatrix withinSubjectContrast(double[] factor1Values, 
			double[] factor2Values)
	throws IllegalArgumentException
	{
		if (factor1Values == null || factor1Values.length < 2 ||
				factor2Values == null || factor2Values.length < 2)
			throw new IllegalArgumentException("must specify at least 2 values for each factor");

		RealMatrix factor1Contrast = withinSubjectContrast(factor1Values);
		RealMatrix factor2Contrast = withinSubjectContrast(factor2Values);
		
		return MatrixUtils.KroneckerProduct(factor1Contrast, factor2Contrast);
	}
	
	/**
	 * Create a within subject contrast (U) for polynomial trends
	 * in three factors using a Kronecker product. It assumes Factor 1
	 * varies most slowly and that Factor 3 varies most rapidly.      
	 * 
	 * @param factor1Values values of first factor
	 * @param factor2Values values of second factor
	 * @param factor3Values values of third factor
	 * @return 3 factor polynomial contrast matrix
	 * @throws IllegalArgumentException
	 */
	public static RealMatrix withinSubjectContrast(double[] factor1Values,
			double[] factor2Values, double[] factor3Values)
	{
		if (factor1Values == null || factor1Values.length < 2 ||
				factor2Values == null || factor2Values.length < 2 || 
				factor3Values == null || factor2Values.length < 3)
			throw new IllegalArgumentException("must specify at least 2 values for each factor");
		
		RealMatrix factor1Contrast = withinSubjectContrast(factor1Values);
		RealMatrix factor2Contrast = withinSubjectContrast(factor2Values);
		RealMatrix factor3Contrast = withinSubjectContrast(factor2Values);
		
		RealMatrix korneckerProduct12 = 
			MatrixUtils.KroneckerProduct(factor1Contrast, factor2Contrast);
		return MatrixUtils.KroneckerProduct(korneckerProduct12, factor3Contrast);
	}
	
}
