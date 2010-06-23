package edu.cudenver.bios.power;

public abstract class Power
{	
	String errMsg = null;
	
	double nominalPower;
	
	double actualPower; 
	
	double totalSampleSize;
	
	double alpha;
	
	public Power(double nominalPower, double actualPower, double totalSampleSize, double alpha)
	{
		this.nominalPower = nominalPower;
		this.actualPower = actualPower;
		this.totalSampleSize = totalSampleSize;
		this.alpha = alpha;
	}
	
	public abstract String toXML();
	
	public double getNominalPower()
	{
		return nominalPower;
	}

	public void setNominalPower(double nominalPower)
	{
		this.nominalPower = nominalPower;
	}

	public double getActualPower()
	{
		return actualPower;
	}

	public void setActualPower(double actualPower)
	{
		this.actualPower = actualPower;
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
