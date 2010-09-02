package edu.cudenver.bios.power.test.general;

import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters.MeanPair;
import junit.framework.TestCase;

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
