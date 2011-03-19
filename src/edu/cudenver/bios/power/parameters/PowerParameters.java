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
//	// wrapper class to allow easy iteration / peek into 
//	// a list
//    protected class PeekableList<T>
//    {
//        ArrayList<T> list = new ArrayList<T>();
//        T currentItem = null;
//        int currentIndex = -1;
//        
//        public void add(T object) { list.add(object); }
//        
//        public T first()
//        {
//            if (list.size() > 0)
//            {
//                currentItem = list.get(0);
//                currentIndex = 0;
//            }
//            return currentItem;
//        }
//        
//        public T next()
//        {
//            if (list.size() > 0 && currentIndex < list.size()-1)
//            {
//                currentIndex++;
//                currentItem = list.get(currentIndex);
//            }
//            else
//            {
//                currentItem = null;
//            }
//            return currentItem;
//        }
//        
//        public T current()
//        {
//            return currentItem;
//        }
//        
//        public int size() { return list.size(); }
//    }
    
	/**
	 * List of power values
	 */
    ArrayList<Double> powerList = new ArrayList<Double>();
    
	/**
	 * List of sample size values
	 */
    ArrayList<Integer> sampleSizeList = new ArrayList<Integer>();
    
    /**
     * power level (type I error level)
     */
    ArrayList<Double> alphaList = new ArrayList<Double>();
    
    /**
     * Create an empty power parameter object
     */
    public PowerParameters() {}

    /**
     * Get the list of desired power values
     * @return list of desired power values
     */
	public ArrayList<Double> getPowerList()
	{
		return powerList;
	}

	/**
	 * Get the list of sample sizes.
	 * @return list of sample sizes
	 */
	public ArrayList<Integer> getSampleSizeList()
	{
		return sampleSizeList;
	}

	/**
	 * Get the list of Type I error (alpha) values.
	 * @return list of alpha values
	 */
	public ArrayList<Double> getAlphaList()
	{
		return alphaList;
	}
    
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

}
