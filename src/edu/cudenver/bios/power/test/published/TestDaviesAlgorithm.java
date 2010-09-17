package edu.cudenver.bios.power.test.published;

import java.util.ArrayList;

import edu.cudenver.bios.power.glmm.WeightedSumOfNoncentralChiSquaresDistribution;
import edu.cudenver.bios.power.glmm.ChiSquareTerm;
import junit.framework.TestCase;

public class TestDaviesAlgorithm extends TestCase
{
	private static final double ACCURACY = 0.001;
	private static final double[] QUANTILES = {2,4,6,8,10};
	private static final double[] lambdas = {7,-3,5};
	private static final double[] omegas = {1,2,1};
	private static final double[] degreesOfFreedom = {10,10,10};
	
	private ArrayList<ChiSquareTerm> chiSquareTerms = new ArrayList<ChiSquareTerm>();
	
	WeightedSumOfNoncentralChiSquaresDistribution dist = null;
	public void setUp()
	{
		// build the chi-square terms
		for(int i = 0; i < lambdas.length && i < omegas.length && i < degreesOfFreedom.length; i++)
		{
			chiSquareTerms.add(new ChiSquareTerm(lambdas[i], degreesOfFreedom[i], omegas[i]));
		}

		dist	= new WeightedSumOfNoncentralChiSquaresDistribution(chiSquareTerms, ACCURACY);
	}
	
	public void testQuantile()
	{
		for(double quantile: QUANTILES)
		{
			double p = dist.cdf(quantile);
			System.out.println(quantile + "\t" + p);
		}
	}
}
