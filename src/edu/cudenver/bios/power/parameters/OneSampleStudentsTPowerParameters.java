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
package edu.cudenver.bios.power.parameters;

import java.util.ArrayList;

/**
 * Input parameters for the one sample students' t test
 * @author Sarah Kreidler
 *
 */
public class OneSampleStudentsTPowerParameters extends PowerParameters
{
	// class to hold null and alternative means
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
    boolean twoTailed = true;
    
    /**
     * Constructor.
     */
    public OneSampleStudentsTPowerParameters()
    {
    	super();
    }
    
    /**
     * Add a standard deviation to the list of power calculations
     * 
     * @param sigma standard deviation
     */
    public void addVariance(double sigma)
    {
    	sigmaList.add(new Double(sigma));    	
    }
    
    /**
     * Add a set of null and alternative means to the list of power calculations
     * @param mu0 estimated mean under the null hypothesis 
     * @param muA estimated mean under the alternative hypothesis
     */
    public void addMeans(double mu0, double muA) 
    {
    	if (powerList.size() > 0 && sampleSizeList.size() > 0) 
    		throw new IllegalArgumentException("Must leave one of power, sample size, or null/aternative mean blank");
    	
    	meansList.add(new MeanPair(mu0, muA));
    }

    /**
     * Indicates if a two tailed calculation is being used
     * @return true if two tailed, false if one tailed
     */
	public boolean isTwoTailed()
	{
		return twoTailed;
	}

	/**
	 * Set whether a two tailed power should be calculated
	 * @param twoTailed if true, uses a two-tailed calculation, otherwise one-tailed
	 */
	public void setTwoTailed(boolean twoTailed)
	{
		this.twoTailed = twoTailed;
	}

	/**
	 * Get the list of mean (null and alternative) pairs
	 * @return list of null and alternative means
	 */
	public ArrayList<MeanPair> getMeansList()
	{
		return meansList;
	}

	/**
	 * Get list of variance values
	 * @return list of variance values
	 */
	public ArrayList<Double> getSigmaList()
	{
		return sigmaList;
	}

	
}
