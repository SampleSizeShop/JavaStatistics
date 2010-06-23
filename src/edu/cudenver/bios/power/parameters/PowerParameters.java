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
    
	public Double getFirstAlpha()
	{
	    return alphaList.first();
	}
	
    public Double getNextAlpha()
    {
        return alphaList.next();
    }
    
    public Double getCurrentAlpha()
    {
        return alphaList.current();
    }
    
    public Integer getFirstSampleSize()
    {
        return sampleSizeList.first();
    }
    
    public Integer getNextSampleSize()
    {
        return sampleSizeList.next();
    }
    
    public Integer getCurrentSampleSize()
    {
        return sampleSizeList.current();
    }
    
    public Double getFirstPower()
    {
        return powerList.first();
    }
    
    public Double getNextPower()
    {
        return powerList.next();
    }
    
    public Double getCurrentPower()
    {
        return powerList.current();
    }
}
