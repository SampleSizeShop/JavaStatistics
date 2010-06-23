package edu.cudenver.bios.power;

public class GLMMPower extends Power
{
	double betaScale;
	double sigmaScale;
	
	public GLMMPower(double alpha, double nominalPower, double actualPower, int sampleSize,
			double betaScale, double sigmaScale)
	{
		super(nominalPower, actualPower, sampleSize, alpha);

		this.betaScale = betaScale;
		this.sigmaScale = sigmaScale;
	}
	
	public double getBetaScale()
	{
		return betaScale;
	}

	public void setBetaScale(double betaScale)
	{
		this.betaScale = betaScale;
	}

	public double getSigmaScale()
	{
		return sigmaScale;
	}

	public void setSigmaScale(double sigmaScale)
	{
		this.sigmaScale = sigmaScale;
	}

	public String toXML()
	{
		return "<power alpha='"+ alpha +"' nominalPower='" + nominalPower + "' actualPower='"+ actualPower +
		"' sampleSize='"+ totalSampleSize + "' betaScale='"+ betaScale +"' sigmaScale='"+ sigmaScale +"' />";
	}
}
