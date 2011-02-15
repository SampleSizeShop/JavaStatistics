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
package edu.cudenver.bios.power.test.published;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test suite for the unit tests generating data for reproducible
 * research purposes
 * 
 * @author Sarah Kreidler
 *
 */
public class GLMMFixedTestSuite
{
    public static Test suite() {

        TestSuite suite = new TestSuite();
  
        suite.addTestSuite(TestConditionalUnivariate.class);
        suite.addTestSuite(TestConditionalMultivariate.class);
        suite.addTestSuite(TestHotellingLawleyApproximateQuantile.class);
        suite.addTestSuite(TestHotellingLawleyApproximateUnconditional.class);
        suite.addTestSuite(TestHotellingLawleyExactQuantile.class);
        suite.addTestSuite(TestHotellingLawleyExactUnconditional.class);
        suite.addTestSuite(TestUnirepApproximateQuantile.class);
        suite.addTestSuite(TestUnirepApproximateUnconditional.class);
        suite.addTestSuite(TestUnirepExactQuantile.class);
        suite.addTestSuite(TestUnirepExactUnconditional.class);

        return suite;
    }

    /**
     * Runs the test suite using the textual runner.
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
