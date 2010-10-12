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

import edu.cudenver.bios.power.glmm.NonCentralFDistribution;
import junit.framework.TestCase;

/**
 * Unit test for non-central F implementation using Tiku approximation
 * @author Sarah Kreidler
 *
 */
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
