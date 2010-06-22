package edu.cudenver.bios.power.parameters;

import java.util.ArrayList;
import java.util.List;

import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;

/**
* Container class for power input parameters.  Most statistical tests
* will need to extend this class.
* 
* @author Sarah Kreidler
*
*/
public abstract class PowerParameters
{
	/**
	 * List of power values
	 */
	ArrayList<Double> powerList = new ArrayList<Double>();

	/**
	 * List of sample size values
	 */
    ArrayList<Integer> sampleSizeList = new ArrayList<Integer>();
    
    /**
     * Alpha level (type I error level)
     */
    ArrayList<Double> alphaList = new ArrayList<Double>();
    
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

	public ArrayList<Double> getPowerList()
	{
		return powerList;
	}

	public ArrayList<Integer> getSampleSizeList()
	{
		return sampleSizeList;
	}

	public ArrayList<Double> getAlphaList()
	{
		return alphaList;
	}
    
    
}
