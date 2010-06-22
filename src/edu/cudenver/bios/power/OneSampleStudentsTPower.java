package edu.cudenver.bios.power;

public class OneSampleStudentsTPower extends Power
{
	double sigma;
	double detectableDifference;
	
	public OneSampleStudentsTPower(double alpha, double power, int sampleSize, 
			double detectableDifference, double sigma)
	{
		super(power, sampleSize, alpha);
		this.sigma = sigma;
		this.detectableDifference = detectableDifference;
	}
	
	public double getSigma()
	{
		return sigma;
	}

	public void setSigma(double sigma)
	{
		this.sigma = sigma;
	}

	public double getDetectableDifference()
	{
		return detectableDifference;
	}

	public void setDetectableDifference(double detectableDifference)
	{
		this.detectableDifference = detectableDifference;
	}

	public String toXML()
	{
		return "<power alpha='"+ alpha +"' power='"+power+"' sampleSize='"+ totalSampleSize +
			"' difference='"+detectableDifference +"' variance='"+ sigma +"'>";
	}
}
