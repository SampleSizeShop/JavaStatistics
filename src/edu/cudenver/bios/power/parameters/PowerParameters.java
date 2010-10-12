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
* Container class for power input parameters.  Most statistical tests
* will need to extend this class.
* 
* @author Sarah Kreidler
*
*/
public abstract class PowerParameters
{
	// wrapper class to allow easy iteration / peek into 
	// a list
    protected class PeekableList<T>
    {
        ArrayList<T> list = new ArrayList<T>();
        T currentItem = null;
        int currentIndex = -1;
        
        public void add(T object) { list.add(object); }
        
        public T first()
        {
            if (list.size() > 0)
            {
                currentItem = list.get(0);
                currentIndex = 0;
            }
            return currentItem;
        }
        
        public T next()
        {
            if (list.size() > 0 && currentIndex < list.size()-1)
            {
                currentIndex++;
                currentItem = list.get(currentIndex);
            }
            else
            {
                currentItem = null;
            }
            return currentItem;
        }
        
        public T current()
        {
            return currentItem;
        }
        
        public int size() { return list.size(); }
    }
    
	/**
	 * List of power values
	 */
    PeekableList<Double> powerList = new PeekableList<Double>();
    
	/**
	 * List of sample size values
	 */
    PeekableList<Integer> sampleSizeList = new PeekableList<Integer>();
    
    /**
     * power level (type I error level)
     */
    PeekableList<Double> alphaList = new PeekableList<Double>();
    
    /**
     * Create an empty power parameter object
     */
    public PowerParameters() {}

    /**
     *  Add a power to the list of calculations
     *  
     *  @param power power to include in the sample size calculations
     */
    public void addPower(double power)
    {
    	powerList.add(new Double(power));
    }
    
    /**
     *  Add a sample size to the list of calculations
     *  
     *  @param sampleSize total sample size to include in the power calculation
     */
    public void addSampleSize(int sampleSize)
    {
    	sampleSizeList.add(new Integer(sampleSize));    	
    }

    /**
     * Add an alpha level to the list of calculations
     * 
     * @param alpha type I error rate for the power or sample size calculation
     */
    public void addAlpha(double alpha)
    {
    	alphaList.add(new Double(alpha));
    }
    
    /**
     * Iterate to the first type I error in the list
     * @return first alpha
     */
	public Double getFirstAlpha()
	{
	    return alphaList.first();
	}
	
	/**
	 * Iterate to the next type I error in the list
	 * @return next alpha
	 */
    public Double getNextAlpha()
    {
        return alphaList.next();
    }
    
    /**
     * Peek at the currently active alpha value
     * @return current alpha
     */
    public Double getCurrentAlpha()
    {
        return alphaList.current();
    }
    
    /**
     * Iterate to the first (per group) sample size in the list
     * @return first sample size
     */
    public Integer getFirstSampleSize()
    {
        return sampleSizeList.first();
    }
    
    /**
     * Iterate to the next (per group) sample size in the list
     * or null if at the end of the list
     * @return next sample size
     */
    public Integer getNextSampleSize()
    {
        return sampleSizeList.next();
    }
    
    /**
     * Peek at the currently active sample size value
     * @return current sample size
     */
    public Integer getCurrentSampleSize()
    {
        return sampleSizeList.current();
    }
    
    /**
     * Iterate to the first power in the list (specified for sample
     * size or detectable difference calculations)
     * @return first power
     */
    public Double getFirstPower()
    {
        return powerList.first();
    }
    
    /**
     * Iterate to the next power in the list (specified for sample
     * size or detectable difference calculations)
     * @return next power
     */
    public Double getNextPower()
    {
        return powerList.next();
    }
    
    /**
     * Peek at the currently active power value
     * @return current power
     */
    public Double getCurrentPower()
    {
        return powerList.current();
    }
}
