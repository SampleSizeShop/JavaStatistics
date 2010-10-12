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
package edu.cudenver.bios.power.test.general;

import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters.MeanPair;
import junit.framework.TestCase;

/**
 * Unit test for internal list iteration in the PowerParameters object
 * @author Sarah Kreidler
 *
 */
public class TestPowerParameters extends TestCase
{
    private class MyParams extends PowerParameters {}
    private double[] alphaList = {0.1, 0.2, 0.3, 0.4};
    private int[] sampleSizeList = {10, 20, 30, 40};
    
    public void testList()
    {
        MyParams params = buildParams();

        int alphaCount = 0;
        for(Double alpha = params.getFirstAlpha(); alpha != null; alpha = params.getNextAlpha())
        {
            System.out.println("alpha = " + alpha);
            alphaCount++;
        }
        assertEquals(alphaCount, alphaList.length);
        
        int sampleSizeCount = 0;
        for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
            sampleSize = params.getNextSampleSize())
        {
            System.out.println("sample size = " + sampleSize);
            sampleSizeCount++;
        }
        assertEquals(sampleSizeCount, sampleSizeList.length);
    }
    
    public void testNestedList()
    {
        MyParams params = buildParams();

        int count = 0;
        for(Double alpha = params.getFirstAlpha(); alpha != null; alpha = params.getNextAlpha())
        {
            for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
                sampleSize = params.getNextSampleSize())
            {
                System.out.println("alpha = " + alpha + ", sampleSize = " + sampleSize);
                count++;
            }
        }
        assertEquals(count, alphaList.length*sampleSizeList.length);
    }
    
    private MyParams buildParams()
    {
        MyParams params = new MyParams();
        
        for(double alpha: alphaList) params.addAlpha(alpha);
        for(int sampleSize: sampleSizeList) params.addSampleSize(sampleSize);
        
        return params;
    }
}
