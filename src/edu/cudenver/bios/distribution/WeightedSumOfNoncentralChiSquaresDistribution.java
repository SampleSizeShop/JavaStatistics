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
package edu.cudenver.bios.distribution;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;


/**
 * Distribution of a weighted sum of noncentral chi square variables
 * using Davie's algorithm,"Distribution of a linear combination of 
 * non-central chi-squared random variables" (1980), Applied Statistics
 * 
 * @author Sarah Kreidler
 *
 */
public class WeightedSumOfNoncentralChiSquaresDistribution 
{
	// lambda coefficient on the normal term added onto the sum of chi-squares
	private static final double DEFAULT_GAUSSIAN_COEFFICIENT = 0.0;
	// default tau-sqaured term
	private static final double DEFAULT_TAU_SQUARED = 0.0;

	// maximum number of iterations allowed
	private static final int MAX_STEPS = 200000;
	// log(2) / 8 - no idea where this comes from, ask Davies
	private static final double LN_2_DIV_8 = Math.log(2.0) / 8.0;
	
	// parameters describing the weights (lambda), degrees of freedom (nu), 
	// and non-centralities (omega) of the chi square variables
	protected List<ChiSquareTerm> chiSquareTerms;
	// map containing the ranks of the chi square terms (ranked by absolute value of lambda)
	HashMap<Double, Integer> chiSquareRankMap = new HashMap<Double, Integer>();
	
	// min and max lambda values
	protected double maxLambda;
	protected double maxLambdaAbsValue;
	protected double minLambda;
	
	// maximum error allowed for the computed probability
	protected double accuracy;
	
	// coefficient of the normal term
	protected double normalCoefficient = DEFAULT_GAUSSIAN_COEFFICIENT;
	
	// simple counter class to maintain threadsafety for the iteration counter
	private class Counter {
		int count;
		int maxSteps;
		public Counter(int maxSteps) 
		{ 
			count = 0; 
			this.maxSteps = maxSteps;
		}
		public void increment() 
		{ 
			count++; 
			if (count > maxSteps)
				throw new RuntimeException("Exceeded max iterations");
		}
		public int getCount() { return count; }
	};
	
	// sorting function for the lambdas
	private class MinMaxComparator implements Comparator<ChiSquareTerm>
	{
		public MinMaxComparator() {};
		@Override
		public int compare(ChiSquareTerm cs1, ChiSquareTerm cs2)
		{
			if (cs1.getLambda() < cs2.getLambda())
				return -1;
			else if (cs1.getLambda() > cs2.getLambda())
				return 1;
			else
				return 0;
		}
	}
	
	// rank sort for the lambdas
	private class RankOrderComparator implements Comparator<Double>
	{
		@Override
		public int compare(Double cs1, Double cs2)
		{
			double absLambda1 = Math.abs(cs1);
			double absLambda2 = Math.abs(cs2);

			if (absLambda1 < absLambda2)
				return -1;
			else if (absLambda1 > absLambda2)
				return 1;
			else
				return 0;
		}
	}
	
	// container class for probability information
	private class TailProbabilityBound 
	{
		public double bound;
		public double cutoff;
		
		public TailProbabilityBound(double bound, double cutoff)
		{
			this.bound = bound;
			this.cutoff = cutoff;
		}
	}
	
	/**
	 * Create a distribution for the specified weighted sum of non-central chi squares
	 * 
	 * @param chiSquareTerms list of weighted chi-square terms
	 * @param normalCoefficient coefficient on the normal term
	 * @param accuracy maximum error allowed in probability
	 * @throws IllegalArgumentException
	 */
	public WeightedSumOfNoncentralChiSquaresDistribution(List<ChiSquareTerm> chiSquareTerms, 
			double normalCoefficient, double accuracy) 
	throws IllegalArgumentException
	{
		if (chiSquareTerms == null || chiSquareTerms.isEmpty())
			throw new IllegalArgumentException("No chi-square terms specified");

		if (Double.isNaN(normalCoefficient))
			throw new IllegalArgumentException("Invalid coefficient for the normal term");
		
		if (Double.isNaN(accuracy) || accuracy <= 0)
			throw new IllegalArgumentException("Accuracy must be greater than 0");
		
		this.chiSquareTerms = chiSquareTerms;
		this.accuracy = accuracy;
		this.normalCoefficient = normalCoefficient;
		
		// find the min/max lambda  (truncate min at 0)
		maxLambda = (Collections.max(chiSquareTerms, new MinMaxComparator())).getLambda();
		minLambda = (Collections.min(chiSquareTerms, new MinMaxComparator())).getLambda();
		maxLambdaAbsValue = (maxLambda < -minLambda ? -minLambda : maxLambda);
		
		if (maxLambda == 0 && minLambda == 0 && normalCoefficient == 0)
			throw new IllegalArgumentException("At least one of min/max lambda values or coefficient of the normal term must be non-zero");
		
		// Build a list of ranks for chi squares based on the absolute value of the lambdas
		ArrayList<Double> sortedList = new ArrayList<Double>();
		for(ChiSquareTerm chiSquare: this.chiSquareTerms) sortedList.add(chiSquare.getLambda());
		Collections.sort(sortedList, new RankOrderComparator());
		int rank = 0;
		for(Double val: sortedList)
		{
			chiSquareRankMap.put(val, rank);
			rank++;
		}	
	}
	
	/**
	 * Create a distribution for the specified weighted sum of non-central chi squares
	 * assuming a 0 as the coefficient on the normal term
	 * 
	 * @param chiSquareTerms list of weighted chi-square terms
	 * @param accuracy maximum error allowed in probability
	 * @throws IllegalArgumentException
	 */
	public WeightedSumOfNoncentralChiSquaresDistribution(List<ChiSquareTerm> chiSquareTerms, 
			double accuracy) throws IllegalArgumentException
	{
		this(chiSquareTerms, DEFAULT_GAUSSIAN_COEFFICIENT, accuracy); 
	}
	
	/**
	 * Returns the cdf of the specified weighted sum of chi-squares distribution 
	 * at the given quantile.
	 * 
	 * @param quantile point at which the cdf is evaluated: Pr(Q &lt;= quantile)
	 * @return probability Pr(Q &lt;= quantile)
	 */
	public double cdf(double quantile)
	throws RuntimeException
	{						
		double prob = 0;
		// expected value of ???
		double mean = 0;
		// initialize sigma
		double sigmaSquared = normalCoefficient * normalCoefficient;
		// initialize the standard deviation of the normal term
		double sd = sigmaSquared;
		// convergence factor
		double tauSquared;
		
		// 
		// sum over something TODO: what?
		for(ChiSquareTerm chiSquare: chiSquareTerms)
		{
			double lambda = chiSquare.getLambda();
			double df = chiSquare.getDegreesOfFreedom();
			double nonCentrality = chiSquare.getNoncentrality();
			
			if (!(df < 0 || nonCentrality < 0))
			{
				sd += lambda * lambda * (2 * df + 4 * nonCentrality);
				mean += lambda * (df + nonCentrality);
			}
			else
			{
				// TODO: throw runtime exception??
			}
		}

		// if sd term is 0, all probability density is piled on 0.
		// thus, we return prob 1 for quantile=0, 0 otherwise
		if (sd == 0)
		{
			if (quantile != 0)
				return 0;
			else
				return 1;
		}
		sd = Math.sqrt(sd);

		// initialize an interation counter  
		Counter counter = new Counter(MAX_STEPS);

		double halfAccuracy = 0.5 * accuracy;

		double U = findTruncationPoint(16 / sd, sigmaSquared, halfAccuracy, counter);
		if (quantile != 0 || maxLambdaAbsValue > 0.07 * sd)
		{
			// check if a covergence factor is helpful
			double convergenceFactor = calculateConvergenceFactor(quantile, counter);
			tauSquared = 0.25*accuracy/convergenceFactor; 
			if (calculateIntegrationError(U, sigmaSquared, tauSquared, counter) < 0.2 * accuracy)
			{
				sigmaSquared += tauSquared;
				U = findTruncationPoint(16 / sd, sigmaSquared, halfAccuracy, counter);
			}
		}



		System.out.println("sd: " + sd);

		System.out.println("mean: " + mean);
		System.out.println("U: " + U);
		
		// Auxiliary integration loop
		double numTermsMain = 0;
		double numTermsAux = 0;
		double integralSum = 0;
		double integrationLimit = 0;
		double integrationInterval = 0;
		do
		{
			// find the range of the distribution, determine if the quantile falls outside 
			// of the range
			double cutoff = findCutoffPoint(4.5 / sd, mean, sigmaSquared, halfAccuracy, counter);
			double cutoffDiffUpper = cutoff - quantile;
			System.out.println("cutoffDiffUpper: " + cutoffDiffUpper);
			if (cutoffDiffUpper < 0) return 1; // requested quantile past range
			// get the lower cutoff
			cutoff = findCutoffPoint(-4.5 / sd, mean, sigmaSquared, halfAccuracy, counter);
			double cutoffDiffLower = quantile - cutoff;
			System.out.println("cutoffDiffLower: " + cutoffDiffLower);
			if (cutoffDiffLower < 0) return 0;

			// pick the larger potential integration interval
			integrationInterval = ((cutoffDiffUpper > cutoffDiffLower) ? cutoffDiffUpper : cutoffDiffLower);
			integrationInterval = 2 * Math.PI / integrationInterval;

			// calculate the #terms required for main and auxilliary integrations
			numTermsMain = U/integrationInterval;
			numTermsAux = 3.0/Math.sqrt(halfAccuracy);
			integralSum = 0;
			if (numTermsMain > 1.5 * numTermsAux)
			{
				// perform the auxilliary integration, provided we have enough iterations left
				if (numTermsAux > MAX_STEPS-counter.getCount())
					throw new RuntimeException("Number of auxiliary integration terms exceeds max number of iteration steps allowed");  

				double integrationIntervalAux = U / numTermsAux;
				integrationLimit = 2*Math.PI/integrationIntervalAux;
				System.out.println("integrationLimit: " + integrationLimit);
				if (integrationLimit <= Math.abs(quantile))
				{
					double lowerConvergenceFactor = calculateConvergenceFactor(quantile-integrationLimit, counter);
					double upperConvergenceFactor = calculateConvergenceFactor(quantile+integrationLimit, counter);

					tauSquared = lowerConvergenceFactor + upperConvergenceFactor;
					tauSquared = (accuracy / 3) / (1.1 * tauSquared);

					accuracy *= 0.67;
					integralSum += integrate((int) Math.round(numTermsAux), integrationIntervalAux, quantile, tauSquared);
					sigmaSquared += tauSquared;

					U = findTruncationPoint(U, sigmaSquared, 0.25*accuracy, counter);
					accuracy *= 0.75;
				}
			}
		}
		while (numTermsMain > 1.5 * numTermsAux && integrationLimit <= Math.abs(quantile) &&
                !Thread.currentThread().isInterrupted());
		
		// perform main integration
		if (numTermsMain > MAX_STEPS-counter.getCount())
			throw new RuntimeException("Number of main integration terms exceeds max number of iteration steps allowed");  
		integralSum += integrate((int) Math.round(numTermsMain), integrationInterval, quantile, Double.NaN);

		prob = 0.5 - integralSum;
		return prob;
	}
	
	/*
	 * Integrate
	 */
	private double integrate(int numTerms, double integrationInterval, double quantile, double TauSquared)
	{
		double value = 0;
		for(int k = numTerms; k >=0; k--)
		{
			double U = (k + 0.5)*integrationInterval;
			double sum1 = -2*U*quantile;
			double sum2 = Math.abs(sum1);
			double sum3 = -0.5*normalCoefficient*normalCoefficient*U*U;
			
			for(ChiSquareTerm chiSquare: chiSquareTerms)
			{
				double lambda = chiSquare.getLambda();
				double df = chiSquare.getDegreesOfFreedom();
				double nonCentrality = chiSquare.getNoncentrality();
				
				double X = 2 * lambda * U;
				double Y = X*X;
				sum3 += -0.25*df*computeLog(Y, true);
				Y = nonCentrality*X / (1+Y);
				double Z = df * Math.atan(X) + Y;
				sum1 += Z;
				sum2 += Math.abs(Z);
				sum3 += -0.5*X*Y;
			}
			
			double partialValue = (integrationInterval/Math.PI)*Math.exp(sum3) / U;
			if (!Double.isNaN(TauSquared))
				partialValue *= (1 - Math.exp(-0.5*TauSquared*U*U));
			
			sum1 = Math.sin(0.5*sum1)*partialValue;
			sum2 = 0.5*sum2*partialValue;
			// TODO: return sum1 from auxiliary integration
			
			value += sum1;
		}
		
		return value;
	}
	
	/**
	 * Find the upper cutoff point for the integral such that P(Q > c) &lt; accuracy
	 * for U > 0, and 
	 * @param startCutoff
	 * @param mean
	 * @param sigmaSquared
     * @param acc
	 * @param counter
	 * @return
	 */
	private double findCutoffPoint(double startCutoff, double mean, double sigmaSquared, double acc, Counter counter)
	throws RuntimeException
	{
		double u2 = startCutoff;
		double u1 = 0;
		double c1 = mean;
		double c2;
		double rb = minLambda;
		
		if (u2 > 0) rb = maxLambda;
		
		rb *= 2;
		
		double u = u2 / (1 + u2*rb);
		TailProbabilityBound bound = findTailProbabilityBound(u, sigmaSquared, counter);
		c2 = bound.cutoff;
		while (bound.bound > acc && !Thread.currentThread().isInterrupted())
		{
			u1 = u2;
			c2 = bound.cutoff;
			c1 = c2;
			u2 *=2;
			u = u2 / (1 + u2*rb);
			bound = findTailProbabilityBound(u, sigmaSquared, counter);
		}
		u = (c1 - mean)/(bound.cutoff - mean);

		while (u < 0.9 && !Thread.currentThread().isInterrupted())
		{
			u = (u1 + u2)/2;
			bound = findTailProbabilityBound(u/(1+ u*rb), sigmaSquared, counter);
			if (bound.bound > acc)
			{
				u1 = u;
				c1 = bound.cutoff;
			}
			else
			{
				u2 = u;
				c2 = bound.cutoff;
			}
			u = (c1 - mean)/(c2 - mean);
		}
		return c2;
	}
	
	/*
	 * Bound on the tail probability for a specified cutoff
	 */
	private TailProbabilityBound findTailProbabilityBound(double startU, double sigmaSquared, Counter counter)
	throws RuntimeException
	{
		counter.increment();
		double U = startU;
		double cutoff = U*sigmaSquared;
		double sum = U*cutoff;
		U *= 2;
		
		for(ChiSquareTerm chiSquare: chiSquareTerms)
		{
			double lambda = chiSquare.getLambda();
			double df = chiSquare.getDegreesOfFreedom();
			double nonCentrality = chiSquare.getNoncentrality();
			
			double X = U * lambda;
			double Y = 1 - X;
			cutoff += lambda*(nonCentrality/Y+df) / Y;
			sum += nonCentrality * ((X/Y)*(X/Y)) + df *( (X*X)/Y + computeLog(-X, false));
		}
		
		TailProbabilityBound bound = new TailProbabilityBound(Math.exp(-0.5 * sum), cutoff);
		
		return bound;
	}
	
	/**
	 * Find truncation point (U) such that the integration error is within the 
	 * user specified accuracy, and that the integration error of U/1.2 is 
	 * greater than the accuracy
	 * 
	 * @param startU starting value for the truncation point
	 */
	private double findTruncationPoint(double startU, double sigmaSquared, double acc, Counter counter)
	throws RuntimeException
	{
		double U = startU / 4;
		double[]  divisors = {2.0, 1.4, 1.2, 1.1}; 
		// find the minimum U such that integration error will be below
		// the specified accuracy
		if (calculateIntegrationError(U, sigmaSquared, DEFAULT_TAU_SQUARED, counter) <= acc)
		{
			// our first guess was too high, so decrease U until we hit 
			// the minimum value with appropriate accuracy
			double Utemp = U;
			U /= 4;
			while(calculateIntegrationError(U, sigmaSquared, DEFAULT_TAU_SQUARED, counter) <= acc &&
                    !Thread.currentThread().isInterrupted())
			{
				Utemp = U;
				U /= 4;
			}
			U = Utemp;
		}
		else
		{
			// our first guess was too low, so increase U until we get
			// good enough accuracy
			do
			{
				U *= 4;
			} 
			while(calculateIntegrationError(U, sigmaSquared, DEFAULT_TAU_SQUARED, counter) > acc &&
                    !Thread.currentThread().isInterrupted());
		}
		
		// now ensure that the truncating error of U/1.2 > accuracy
		double Utemp = U;
		for(double divisor: divisors)
		{
			Utemp /= divisor;
			if (calculateIntegrationError(Utemp, sigmaSquared, DEFAULT_TAU_SQUARED, counter) > acc)
				break;
			U = Utemp;
		}
		
		return U;
	}
	
	/**
	 * Calculate the tau convergence factor
	 * @return convergence factor
	 */
	private double calculateConvergenceFactor(double quantile, Counter counter)
	throws RuntimeException
	{
		counter.increment();

		// sign of the quantile (1, 0, or -1)
		double quantileSign = 0;
		if (quantile > 0)
			quantileSign = 1;
		else if (quantile < 0)
		{
			quantileSign = -1;
		}
		
		double sum = 0;
		double axl = Math.abs(quantile);

		for(int j = chiSquareTerms.size()-1; j >= 0; j--)
		{
			ChiSquareTerm chiSquare = chiSquareTerms.get(j);
			int rankIndex = this.chiSquareRankMap.get(chiSquare.getLambda());
			ChiSquareTerm rankChiSquare = chiSquareTerms.get(rankIndex);
		
			double lambda = rankChiSquare.getLambda();
			double lambdaAbsoluteValue = Math.abs(lambda);
			double df = rankChiSquare.getDegreesOfFreedom();
			double nonCentrality = rankChiSquare.getNoncentrality();
			
			if (lambda*quantileSign > 0)
			{
				double axl1 = axl - lambdaAbsoluteValue*(df + nonCentrality);
				double axl2 = lambda / LN_2_DIV_8;
				
				if (axl1 > axl2)
				{
					axl = axl1;
				}
				else
				{
					if (axl > axl2) axl = axl2;
					
					sum = (axl - axl1)/lambda;
					
					for(int k = j-1; k >=0; k--)
					{
						ChiSquareTerm chiSquareK = chiSquareTerms.get(k);
						int rankIndexK = this.chiSquareRankMap.get(chiSquareK.getLambda());
						ChiSquareTerm rankChiSquareK = chiSquareTerms.get(rankIndexK);

						sum += rankChiSquareK.getDegreesOfFreedom() + rankChiSquareK.getNoncentrality();
					}
					break;
				}
			}
		}
		
//		if (sum > 100)
//			throw new RuntimeException("Summation to big, bad roundoff error"); // TODO: do we need this check, another precision thing?
//		
		return Math.pow(2, sum/4) / (Math.PI*axl*axl);
	}
	
	/**
	 * Calculates the integration error associated with the given
	 * truncation point, U.
	 * 
	 * @return integration error
	 */
	private double calculateIntegrationError(double U, double sigmaSquared, double tauSquared, Counter counter)
	throws RuntimeException
	{
		counter.increment();
		
		double sum1 = 0;
		double sum2 = (sigmaSquared + tauSquared) * U * U;
		double prod1 = 2 * sum2;
		double prod2 = 0;
		double prod3 = 0;
		double ns = 0;
		double X = 0;
		U *= 2;

		for(ChiSquareTerm chiSquare: chiSquareTerms)
		{
			double lambda = chiSquare.getLambda();
			double df = chiSquare.getDegreesOfFreedom();
			double nonCentrality = chiSquare.getNoncentrality();
			
			X = (U*lambda) * (U*lambda);
			sum1 += nonCentrality * X / (1+X);
			if (X <= 1)
			{
				prod1 += df * computeLog(X, true);
			}
			else
			{
				prod2 += df * Math.log(X);
				prod3 += df * computeLog(X, true);
				ns += df;
			}
		}
		
		sum1 *= 0.5;
		prod2 += prod1;
		prod3 += prod1;
		X = Math.exp(-sum1 - 0.25*prod2) / Math.PI;
		double Y = Math.exp(-sum1 - 0.25*prod3) / Math.PI;
		
		// now that we've accumulated terms, build the three forms of
		// the truncation error
		double error1;
		double error2;
		
		if (ns != 0)
		{
			error1 = X * 2 / ns;
		}
		else
		{
			error1 = 1;
		}
		
		if (prod3 > 1)
		{
			error2 = 2.5 * Y;
		}
		else
		{
			error2 = 1;
		}
		
		// take the smaller of the two errors
		if (error2 < error1) error1 = error2;
		
		X = 0.5*sum2;
		
		if (X <= Y)
			error2 = 1;
		else
			error2 = Y/X;
		
		return (error1 < error2 ? error1 : error2);
	}
	
	
	/**
	 * convenience function for log
	 */
	private double computeLog(double x, boolean first)
	{
		if (Math.abs(x) > 0.1)
		{
			return (first ? Math.log(1 + x) : (Math.log(1 + x) - x));
		}
		else
		{
			// deal with values very close to 0
			double y = x / (2 + x);
			double term = 2 * y * y * y; // 2 * Y^3
			double ak = 3.0;
			double s = 2;
			if (!first) s = -x;
			
			s *= y;
			y *= y;
			double s1 = s + term/ak;
			
			while(s1 != s && !Thread.currentThread().isInterrupted())
			{
				ak = ak +2;
				term *= y;
				s = s1;
				s1 = s + term/ak;
			}
			return s;
		}
	}

}
