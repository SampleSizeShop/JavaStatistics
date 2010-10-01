package edu.cudenver.bios.power.test.published;

import junit.framework.Test;
import junit.framework.TestSuite;

public class PublishableResultsTestSuite
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
