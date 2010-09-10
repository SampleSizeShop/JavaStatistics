package edu.cudenver.bios.power.glmm;

import java.util.Collections;
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
	private static final double LN_2_DIV_8 = Math.log10(2.0) / 8.0;
	
	// parameters describing the weights (lambda), degrees of freedom (nu), 
	// and non-centralities (omega) of the chi square variables
	protected List<Double> lambdaList;
	protected List<Double> degreesOfFreedomList;
	protected List<Double> noncentralityList;

	// min and max lambda values
	protected double maxLambda;
	protected double minLambda;
	
	// maximum error allowed for the computed probability
	protected double accuracy;
	
	// coefficient of the normal term
	protected double normalCoefficient = DEFAULT_GAUSSIAN_COEFFICIENT;
	
	// simple counter class to maintain threadsafety for the iteration counter
	private class Counter {
		int count;
		public Counter() { count = 0; }
		public void increment() { count++; }
		public int getCount() { return count; }
	};
	
	/**
	 * Create a distribution for the specified weighted sum of non-central chi squares
	 * 
	 * @param lambda Nx1 vector of coefficients
	 * @param nu Nx1 vector of degrees of freedom
	 * @param omega Nx1 vector of noncentralities
	 * @param normalCoefficient coefficient on the normal term
	 * @param accuracy maximum error allowed in probability
	 * @throws IllegalArgumentException
	 */
	public WeightedSumOfNoncentralChiSquaresDistribution(List<Double> lambda, 
			List<Double> nu, List<Double> omega, double normalCoefficient, double accuracy) 
	throws IllegalArgumentException
	{
		double N = lambda.size();
		
		// make sure all three lists are the same size
		if (nu.size() != N || omega.size() != N)
			throw new IllegalArgumentException("length of lambda, nu, and omega inputs must be equal to one another");
		
		this.lambdaList = lambda;
		this.degreesOfFreedomList = nu;
		this.noncentralityList = omega;
		this.accuracy = accuracy;
		this.normalCoefficient = normalCoefficient;
		
		maxLambda = (Double) Collections.max(lambdaList);
		minLambda = (Double) Collections.min(lambdaList);
		if (maxLambda < -minLambda) maxLambda = -minLambda;
		
		if (maxLambda == 0 && minLambda == 0 && normalCoefficient == 0)
			throw new IllegalArgumentException("At least one of min/max lambda values or coefficient of the normal term must be non-zero");
	}
	
	/**
	 * Create a distribution for the specified weighted sum of non-central chi squares
	 * assuming a 0 as the coefficient on the normal term
	 * 
	 * @param lambda Nx1 vector of coefficients
	 * @param nu Nx1 vector of degrees of freedom
	 * @param omega Nx1 vector of noncentralities
	 * @param accuracy maximum error allowed in probability
	 * @throws IllegalArgumentException
	 */
	public WeightedSumOfNoncentralChiSquaresDistribution(List<Double> lambda, 
			List<Double> nu, List<Double> omega, double accuracy) throws IllegalArgumentException
	{
		this(lambda, nu, omega, DEFAULT_GAUSSIAN_COEFFICIENT, accuracy); 
	}
	
	/**
	 * Returns the cdf of the specified weighted sum of chi-squares distribution 
	 * at the given quantile.
	 * 
	 * @param quantile point at which the cdf is evaluated: Pr(Q &lt;= quantile)
	 * @param lambda Nx1 vector of coefficients
	 * @param nu Nx1 vector of degrees of freedom
	 * @param omega Nx1 vector of noncentralities
	 * @param accuracy maximum error allowed in probability
	 * @throws RuntimeException thrown when the probability cannot be computed within the specified accuracy
	 * @return probability Pr(Q &lt;= quantile)
	 */
	public double cdf(double quantile)
	throws RuntimeException
	{
		// Chi-square values must be positive, so probability is 0 
		// for non-positive values TODO: double check this??
		//if (quantile < 0) return 0;
				
		// make sure at least one of min, max lambda and normal coefficient 
		// are non-zero
		if (minLambda == 0 && maxLambda == 0 && normalCoefficient == 0)
			throw new RuntimeException("Max/min lambda and coefficient of normal term cannot all be 0");
		
		double prob = 0;
		// number of chi-square terms in the sum
		double numChiSquareTerms = lambdaList.size(); 
		// expected value of ???
		double mean = 0;
		// initialize the standard deviation of the normal term
		double sd = normalCoefficient * normalCoefficient;
		
		// sum over something TODO: what?
		for(int i = 0; i < numChiSquareTerms; i++)
		{
			double lambda = lambdaList.get(i);
			double df = degreesOfFreedomList.get(i);
			double nonCentrality = noncentralityList.get(i);
			
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
		Counter counter = new Counter();

		double halfAccuracy = 0.5 * accuracy;

		double U = findTruncationPoint(16 / sd, halfAccuracy, counter);
		
		if (quantile == 0 || maxLambda <= 0.07 * sd)
		{
			// TODO: check if a covergence factor is helpful
		}
		
		// find the range of the distribution, determine if the quantile falls outside 
		// of the range
		double cutoff = findCutoffPoint(4.5 / sd, halfAccuracy, mean, counter);
		double cutoffDiffUpper = cutoff - quantile;
		if (cutoffDiffUpper < 0) return 1; // requested quantile past range
		// get the lower cutoff
		cutoff = findCutoffPoint(-4.5 / sd, halfAccuracy, mean, counter);
		double cutoffDiffLower = quantile - cutoff;
		if (cutoffDiffLower < 0) return 0;
		
		// pick the larger potential integration interval
		double integrationInterval = ((cutoffDiffUpper > cutoffDiffLower) ? cutoffDiffUpper : cutoffDiffLower);
		integrationInterval = 2 * Math.PI / integrationInterval;
		
		// calculate the #terms required for main and auxilliary integrations
		double numTermsAux = U/integrationInterval;
		double numTermsMain = 3.0/Math.sqrt(halfAccuracy);

		if (numTermsAux > 1.5 * numTermsMain)
		{
			if (numTermsMain > MAX_STEPS)
				throw new RuntimeException("Number of integration terms exceeds max number of iteration steps allowed");  
			
			integrationInterval = U / numTermsMain;
			numTermsMain = Math.round(numTermsMain);
			//X = 2 * Math.PI / integrationInterval;
		}
		
		// perform main integration
		
		return prob;
	}
	
	
	private double integrate(int numTerms, double integrationInterval, double truncationPoint, 
			double quantile, double TauSquared)
	{
		double value = 0;
		for(int k = numTerms; k >=0; k--)
		{
			double U = (k + 0.5)*integrationInterval;
			double sum1 = -2*U*quantile;
			double sum2 = Math.abs(sum1);
			double sum3 = -0.5*normalCoefficient*normalCoefficient*U*U;
			
			for(int j = lambdaList.size()-1; j >= 0; j--)
			{
				double lambda = lambdaList.get(j);
				double df = degreesOfFreedomList.get(j);
				double omega = noncentralityList.get(j);
				
				double X = 2 * lambda * U;
				double Y = X*X;
				sum3 += -0.25*df*computeLog(Y, true);
				Y = omega*X / (1+Y);
				double Z = df * Math.atan(X) + Y;
				sum1 += Z;
				sum2 += Math.abs(Z);
				sum3 += -0.5*X*Y;
			}
			
			double partialValue = (integrationInterval/Math.PI)*Math.exp(sum3) / U;
			if (Double.isNaN(TauSquared))
				partialValue *= (1 - Math.exp(-0.5*TauSquared*U*U));
			
			sum1 = Math.sin(0.5*sum1)*partialValue;
			sum2 = 0.5*sum2*partialValue;
			// TODO: return sum1 from auxilliary integration
			
			value += sum2;
		}
		
		return value;
	}
	
	/**
	 * Find the upper cutoff point for the integral such that P(Q > c) &lt; accuracy
	 * for U > 0, and 
	 * @param acc
	 * @param mean
	 * @param limit
	 * @param counter
	 * @return
	 */
	private double findCutoffPoint(double startCutoff, double acc, double mean, Counter counter)
	{
		
		
		return 0.0;
	}
	
	/**
	 * Find truncation point (U) such that the integration error is within the 
	 * user specified accuracy, and that the integration error of U/1.2 is 
	 * greater than the accuracy
	 * 
	 * @param startU starting value for the truncation point
	 */
	private double findTruncationPoint(double startU, double acc, Counter counter)
	throws RuntimeException
	{
		double U = startU / 4;
		double[]  divisors = {2.0, 1.4, 1.2, 1.1}; 
		// find the minimum U such that integration error will be below
		// the specified accuracy
		if (findIntegrationError(U, DEFAULT_TAU_SQUARED, counter) <= acc)
		{
			// our first guess was too high, so decrease U until we hit 
			// the minimum value with appropriate accuracy
			double Utemp = U;
			U /= 4;
			while(findIntegrationError(U, DEFAULT_TAU_SQUARED, counter) <= acc) 
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
			U *= 4;
			while(findIntegrationError(U, DEFAULT_TAU_SQUARED, counter) > acc) U *= 4;
		}
		
		// now ensure that the truncating error of U/1.2 > accuracy
		for(double divisor: divisors)
		{
			double Utemp = U / divisor;
			if (findIntegrationError(Utemp, DEFAULT_TAU_SQUARED, counter) > acc)
			{
				U = Utemp;
				break;
			}
		}
		
		return U;
	}
	
	/**
	 * Calculates the integration error associated with the given
	 * truncation point, U.
	 * 
	 * @return integration error
	 */
	private double findIntegrationError(double U, double tauSquared, Counter counter)
	throws RuntimeException
	{
		counter.increment();
		if (counter.getCount() > MAX_STEPS)
			throw new RuntimeException("Exceeded max iterations");
		
		double sum1 = 0;
		double sum2 = (normalCoefficient * normalCoefficient + tauSquared) * U * U;
		double prod1 = 2 * sum2;
		double prod2 = 0;
		double prod3 = 0;
		double ns = 0;
		double X = 0;
		U *= 2;
		double N = lambdaList.size();
		for(int j = 0; j < N; j++)
		{
			double df = degreesOfFreedomList.get(j);
			double lambda = lambdaList.get(j);
			double omega = noncentralityList.get(j);
			
			X = (U*lambda) * (U*lambda);
			sum1 += omega * X / (1+X);
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
		X = exponentiate(-sum1 - 0.25*prod2) / Math.PI;
		double Y = exponentiate(-sum1 - 0.25*prod3) / Math.PI;
		
		// now that we've accumulated terms, build the three forms of
		// the truncation error
		double error1;
		double error2;
		double error3;
		
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
	 * 
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
			
			while(s1 != s)
			{
				ak = ak +2;
				term *= y;
				s = s1;
				s1 = s + term/ak;
			}
			return s;
		}
	}
	
	private double exponentiate(double x)
	{
		if (x < -706)
			return 0; // TODO: is this an algol precision artifact?
		else
			return Math.exp(x);
	}
}
