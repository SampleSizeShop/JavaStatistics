package edu.cudenver.bios.power.glmm;

public class ChiSquareTerm
{
		protected double lambda;
		protected double noncentrality;
		protected double df;
		
		public ChiSquareTerm(double lambda, double df, double noncentrality)
		{
			this.lambda = lambda;
			this.noncentrality = noncentrality;
			this.df = df;
		}
		
		public double getLambda() { return lambda; }
		public double getNoncentrality() { return noncentrality; }
		public double getDegreesOfFreedom() { return df; }
}
