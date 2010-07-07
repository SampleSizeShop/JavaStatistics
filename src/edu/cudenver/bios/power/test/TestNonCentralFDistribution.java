package edu.cudenver.bios.power.test;

import edu.cudenver.bios.power.glmm.NonCentralFDistribution;
import junit.framework.TestCase;

public class TestNonCentralFDistribution extends TestCase
{
    double[] ndfList = {2, 25.3, 3000000};
    double[] ddfList = {5.6, 38.1, 1450000};
    double[] ncList = {0, 5, 3.2};    
    
    public void testCDF()
    {
        int numTests = ndfList.length * ddfList.length * ncList.length;
        int successes = 0;
        for(double ndf: ndfList)
        {
            for(double ddf: ddfList)
            {
                for(double nc: ncList)
                {
                    NonCentralFDistribution dist = new NonCentralFDistribution(ndf, ddf, nc);
                    double prob = dist.cdf(5);
                    System.out.println("Ndf=" + ndf + " ddf=" + ddf + " nc=" + nc + " prob=" + prob);
                    if (prob >= 0 && prob <= 1) successes++;
                }
            }
        }
        
        assertEquals(numTests, successes);
    }
    
    public void testInverseCDF()
    {
        int numTests = ndfList.length * ddfList.length * ncList.length;
        int successes = 0;
        for(double ndf: ndfList)
        {
            for(double ddf: ddfList)
            {
                for(double nc: ncList)
                {
                    NonCentralFDistribution dist = new NonCentralFDistribution(ndf, ddf, nc);
                    double f = dist.inverseCDF(0.5);
                    System.out.println("Ndf=" + ndf + " ddf=" + ddf + " nc=" + nc + " f=" + f);
                    if (f > 0) successes++;
                }
            }
        }
        
        assertEquals(numTests, successes);
    }
}
