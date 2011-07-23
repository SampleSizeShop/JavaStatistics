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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.QRDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;

import edu.cudenver.bios.matrix.OrthogonalPolynomialContrast.ContrastType;
import edu.cudenver.bios.utils.Factor;

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

		RealMatrix z = MatrixUtils.getDiagonalMatrix(qrDecomp.getR());
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
	 * for an arbitrary list of factors.  The returned collection includes
	 * the grand mean contrast, all 1-factor main effect contrasts, and
	 * all possible interaction contrasts
	 * 
	 * @param factorList list of factors, including name and value information
	 * @return polynomial contrast collection.
	 * @throws IllegalArgumentException
	 */
	public static OrthogonalPolynomialContrastCollection withinSubjectContrast(List<Factor> factorList)
	throws IllegalArgumentException
	{		
		if (factorList == null || factorList.size() <= 0)
			throw new IllegalArgumentException("no factors specified");

		ArrayList<RealMatrix> zeroTrendList = new ArrayList<RealMatrix>(factorList.size());
		ArrayList<RealMatrix> factorTrendList = new ArrayList<RealMatrix>(factorList.size());
		for(Factor factor : factorList)
		{
			double[] centered = centerAndScale(factor.getValues());
			RealMatrix poly = OrthogonalPolynomials.orthogonalPolynomialCoefficients(centered,centered.length-1);
			zeroTrendList.add(poly.getColumnMatrix(0));
			factorTrendList.add(poly.getSubMatrix(0, poly.getRowDimension()-1, 1, poly.getColumnDimension()-1));
		}
		
		OrthogonalPolynomialContrastCollection results = new OrthogonalPolynomialContrastCollection();
		/*
		 * We need to create contrasts for every possible combination of
		 * zero and factor trends.  The possible combinations can be  
		 * enumerated by the binary representation of the numbers
		 * 0 - 2^(#factors).  In this case, 0 indicates that the zero-trend contrast should
		 * be included in the Kronecker product for a given factor, and 1 indicates that the 
		 * factor trend contrast should be included in the Kronecker product.
		 */
		ArrayList<RealMatrix> kroneckerList = new ArrayList<RealMatrix>(factorList.size());
		ArrayList<Factor> activeFactorList = new ArrayList<Factor>(factorList.size());
		// build the grand mean
		for(RealMatrix zeroTrend : zeroTrendList) kroneckerList.add(zeroTrend);
		results.addContrast(new OrthogonalPolynomialContrast(MatrixUtils.getKroneckerProduct(kroneckerList)));
		// loop over the remaining contrasts
		int totalContrasts = (int) Math.pow(2.0, (double) factorList.size());
		for(int i = 1; i < totalContrasts; i++)
		{
			kroneckerList.clear();
			activeFactorList.clear();
			int mask = 1;
			for(int factorIdx = 0; factorIdx < factorList.size(); factorIdx++, mask *=2)
			{
				if ((i & mask) != 0)
				{
					kroneckerList.add(factorTrendList.get(factorIdx));
					activeFactorList.add(factorList.get(factorIdx));
				}
				else
				{
					kroneckerList.add(zeroTrendList.get(factorIdx));
				}
			}
			// add the appropriate contrast type
			// note that if "i" is a power of 2 then we have a  main effect contrast, else interaction
			results.addContrast(
				new OrthogonalPolynomialContrast(((i & (i-1)) == 0 ? ContrastType.MAIN_EFFECT : 
					ContrastType.INTERACTION),
						activeFactorList, MatrixUtils.getKroneckerProduct(kroneckerList)));
		}
		
		return results;
	}
	
//	/**
//	 * Create a within subject contrast (U) for polynomial trends
//	 * in two factors using a Kronecker product. It assumes Factor 1
//	 *  varies slowly and that Factor 2 varies rapidly.    
//	 * 
//	 * @param factor1Values values of first factor
//	 * @param factor2Values values of second factor
//	 * @return 2 factor polynomial contrast matrix
//	 * @throws IllegalArgumentException
//	 */
//	public static OrthogonalPolynomialContrastCollection withinSubjectContrast(double[] factor1Values, String factor1Name,
//			double[] factor2Values, String factor2Name)
//	throws IllegalArgumentException
//	{
//		if (factor1Values == null || factor1Values.length < 2 ||
//				factor2Values == null || factor2Values.length < 2)
//			throw new IllegalArgumentException("must specify at least 2 values for each factor");
//
//		double[] centeredFactor1 = centerAndScale(factor1Values);
//		double[] centeredFactor2 = centerAndScale(factor2Values);
//		
//		// get the orthogonal polynomials for the 0 - nth polynomials, where n = length of factor array
//		RealMatrix factor1OrthoPoly = orthogonalPolynomialCoefficients(centeredFactor1,centeredFactor1.length-1);
//		RealMatrix factor2OrthoPoly = orthogonalPolynomialCoefficients(centeredFactor2,centeredFactor2.length-1);
//		
//		// extract the zero-order trends
//		RealMatrix zeroTrendFactor1 = factor1OrthoPoly.getColumnMatrix(0);
//		RealMatrix zeroTrendFactor2 = factor2OrthoPoly.getColumnMatrix(0);
//		// extract the 1-nth trends
//		RealMatrix trendsFactor1 = 
//			factor1OrthoPoly.getSubMatrix(0, factor1OrthoPoly.getRowDimension()-1, 
//				1, factor1OrthoPoly.getColumnDimension()-1);
//		RealMatrix trendsFactor2 = 
//			factor2OrthoPoly.getSubMatrix(0, factor2OrthoPoly.getRowDimension()-1, 
//				1, factor2OrthoPoly.getColumnDimension()-1);
//		
//		// build the grand mean
//		RealMatrix grandMean = MatrixUtils.getKroneckerProduct(zeroTrendFactor1, zeroTrendFactor2);
//		
//		// build the main effects contrasts
//		RealMatrix mainEffectFactor1 = MatrixUtils.getKroneckerProduct(trendsFactor1, zeroTrendFactor2);
//		RealMatrix mainEffectFactor2 = MatrixUtils.getKroneckerProduct(zeroTrendFactor1, trendsFactor2);
//		// build the interaction contrast
//		RealMatrix interactionContrast =  MatrixUtils.getKroneckerProduct(trendsFactor1, trendsFactor2);
//		
//		// build the contrast collection
//		OrthogonalPolynomialContrastCollection results = new OrthogonalPolynomialContrastCollection();
//		results.setGrandMean(grandMean);
//		results.addMainEffectContrast(factor1Name, mainEffectFactor1);
//		results.addMainEffectContrast(factor2Name, mainEffectFactor2);
//		//results.addInteractionContrast(factor1Name, factor2Name, interactionContrast);
//		return results;
//	}
//	
//	/**
//	 * Create a within subject contrast (U) for polynomial trends
//	 * in three factors using a Kronecker product. It assumes Factor 1
//	 * varies most slowly and that Factor 3 varies most rapidly.      
//	 * 
//	 * @param factor1Values values of first factor
//	 * @param factor2Values values of second factor
//	 * @param factor3Values values of third factor
//	 * @return 3 factor polynomial contrast matrix
//	 * @throws IllegalArgumentException
//	 */
//	public static OrthogonalPolynomialContrastCollection withinSubjectContrast(double[] factor1Values, String factor1Name,
//			double[] factor2Values, String factor2Name, double[] factor3Values, String factor3Name)
//	{
//		if (factor1Values == null || factor1Values.length < 2 ||
//				factor2Values == null || factor2Values.length < 2 || 
//				factor3Values == null || factor3Values.length < 2)
//			throw new IllegalArgumentException("must specify at least 2 values for each factor");
//		
//		// center and scale the values
//		double[] centeredFactor1 = centerAndScale(factor1Values);
//		double[] centeredFactor2 = centerAndScale(factor2Values);
//		double[] centeredFactor3 = centerAndScale(factor3Values);
//		
//		// get the orthogonal polynomials for the 0 - nth polynomials, where n = length of factor array
//		RealMatrix factor1OrthoPoly = orthogonalPolynomialCoefficients(centeredFactor1,centeredFactor1.length-1);
//		RealMatrix factor2OrthoPoly = orthogonalPolynomialCoefficients(centeredFactor2,centeredFactor2.length-1);
//		RealMatrix factor3OrthoPoly = orthogonalPolynomialCoefficients(centeredFactor3,centeredFactor3.length-1);
//
//		// extract the zero-order trends
//		RealMatrix zeroTrendFactor1 = factor1OrthoPoly.getColumnMatrix(0);
//		RealMatrix zeroTrendFactor2 = factor2OrthoPoly.getColumnMatrix(0);
//		RealMatrix zeroTrendFactor3 = factor3OrthoPoly.getColumnMatrix(0);
//		// extract the 1-nth trends
//		RealMatrix trendsFactor1 = 
//			factor1OrthoPoly.getSubMatrix(0, factor1OrthoPoly.getRowDimension()-1, 
//				1, factor1OrthoPoly.getColumnDimension()-1);
//		RealMatrix trendsFactor2 = 
//			factor2OrthoPoly.getSubMatrix(0, factor2OrthoPoly.getRowDimension()-1, 
//				1, factor2OrthoPoly.getColumnDimension()-1);
//		RealMatrix trendsFactor3 = 
//			factor3OrthoPoly.getSubMatrix(0, factor3OrthoPoly.getRowDimension()-1, 
//				1, factor3OrthoPoly.getColumnDimension()-1);
//		
//		// build the grand mean
//		RealMatrix grandMean = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(zeroTrendFactor1, zeroTrendFactor2),
//					zeroTrendFactor3);
//		
//		// build the main effects contrasts
//		RealMatrix mainEffectFactor1 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(trendsFactor1, zeroTrendFactor2),
//					zeroTrendFactor3);
//		RealMatrix mainEffectFactor2 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(zeroTrendFactor1, trendsFactor2),
//				zeroTrendFactor3);
//		RealMatrix mainEffectFactor3 =  
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(zeroTrendFactor1, zeroTrendFactor2),
//					trendsFactor3);
//		
//		// build the pairwise interaction contrasts
//		RealMatrix interaction12 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(trendsFactor1, trendsFactor2),
//					zeroTrendFactor3);
//		RealMatrix interaction13 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(trendsFactor1, zeroTrendFactor2),
//				trendsFactor3);
//		RealMatrix interaction23 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(zeroTrendFactor1, trendsFactor2),
//					trendsFactor3);
//		
//		// build 3-factor interaction
//		RealMatrix interaction123 = 
//			MatrixUtils.getKroneckerProduct(MatrixUtils.getKroneckerProduct(trendsFactor1, trendsFactor2),
//					trendsFactor3);
//		
//		// build the contrast collection
//		OrthogonalPolynomialContrastCollection results = new OrthogonalPolynomialContrastCollection();
//		results.setGrandMean(grandMean);
//		results.addMainEffectContrast(factor1Name, mainEffectFactor1);
//		results.addMainEffectContrast(factor2Name, mainEffectFactor2);
//		results.addMainEffectContrast(factor3Name, mainEffectFactor3);
////		results.addTwoFactorInteractionContrast(factor1Name, factor2Name, interaction12);
////		results.addTwoFactorInteractionContrast(factor1Name, factor3Name, interaction13);
////		results.addTwoFactorInteractionContrast(factor2Name, factor3Name, interaction23);
////		results.addThreeFactorInteractionContrast(factor1Name, factor2Name, factor3Name, interaction123);
//		return results;
//	}
	
	/**
	 * Center and scale the incoming factor values
	 * @param values
	 */
	private static double[] centerAndScale(double[] values)
	{
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
		return centered;
	}

}
