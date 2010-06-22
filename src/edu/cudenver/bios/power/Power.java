package edu.cudenver.bios.power;

public abstract class Power
{
	double power;
	
	double totalSampleSize;
	
	double alpha;
	
	public Power(double power, double totalSampleSize, double alpha)
	{
		this.power = power;
		this.totalSampleSize = totalSampleSize;
		this.alpha = alpha;
	}
	
	public abstract String toXML();

	public double getPower()
	{
		return power;
	}

	public void setPower(double power)
	{
		this.power = power;
	}

	public double getTotalSampleSize()
	{
		return totalSampleSize;
	}

	public void setTotalSampleSize(double totalSampleSize)
	{
		this.totalSampleSize = totalSampleSize;
	}

	public double getAlpha()
	{
		return alpha;
	}

	public void setAlpha(double alpha)
	{
		this.alpha = alpha;
	}
	
	
}
