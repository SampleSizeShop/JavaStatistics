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
package edu.cudenver.bios.power.test.paper;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test suite for the unit tests generating data for reproducible
 * research purposes
 * 
 * @author Sarah Kreidler
 *
 */
public class PaperTable2TestSuite
{
    public static Test suite() {

        TestSuite suite = new TestSuite();
        // GLMM(F) designs
        suite.addTestSuite(TestConditionalTwoSampleTTest.class);
        suite.addTestSuite(TestConditionalPairedTTest.class);
        suite.addTestSuite(TestConditionalTwoSampleTTest3DPlot.class);
        suite.addTestSuite(TestConditionalUnivariateWithConfidenceLimits.class);
        suite.addTestSuite(TestConditionalMultivariateInteraction.class);
        suite.addTestSuite(TestConditionalMultivariateWithConfidenceLimits.class);
        suite.addTestSuite(TestConditionalOrthogonalPolynomial1Factor.class);
        suite.addTestSuite(TestConditionalOrthogonalPolynomial3Factor.class);
        suite.addTestSuite(TestConditionalOrthogonalPolynomial2Factor.class);

        return suite;
    }

    /**
     * Runs the test suite using the textual runner.
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
