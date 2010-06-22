package edu.cudenver.bios.power.parameters;

import java.util.ArrayList;

import jsc.distributions.StudentsT;

public class OneSampleStudentsTPowerParameters extends PowerParameters
{
	
	public class MeanPair 
	{
		public double mu0 = Double.NaN;
		public double muA = Double.NaN;
		
		public MeanPair(double mu0, double muA)
		{
			this.mu0 = mu0;
			this.muA = muA;
		}
	}
	
    // means under the null and alternative hypotheses
	ArrayList<MeanPair> meansList = new ArrayList<MeanPair>();
    
    // estimated population std dev
    ArrayList<Double> sigmaList = new ArrayList<Double>();
    
    // indicates if a one or two tailed test should be performed
    boolean oneTailed;
    
    public OneSampleStudentsTPowerParameters()
    {
    	super();
    }
    
    /**
     * Add a variance
     */
    public void addVariance(double variance)
    {
    	sigmaList.add(new Double(variance));    	
    }
    
    public void addMeans(double mu0, double muA) 
    {
    	if (powerList.size() > 0 && sampleSizeList.size() > 0) 
    		throw new IllegalArgumentException("Must leave one of power, sample size, or null/aternative mean blank");
    	
    	meansList.add(new MeanPair(mu0, muA));
    }

	public boolean isOneTailed()
	{
		return oneTailed;
	}

	public void setOneTailed(boolean oneTailed)
	{
		this.oneTailed = oneTailed;
	}

	public ArrayList<MeanPair> getMeansList()
	{
		return meansList;
	}

	public ArrayList<Double> getSigmaList()
	{
		return sigmaList;
	}
    
	
}
