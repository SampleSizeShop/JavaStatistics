/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2015 Regents of the University of Colorado.
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
package edu.ucdenver.bios.power;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * Test case for approximate unconditional power for the HLT.  Values should match
 * approximate unconditional power values from Table II in Glueck & Muller 2003.
 * 
 * @author Sarah Kreidler
 *
 */
public class Homework1Cluster04 {

	private static final String DATA_FILE =  "Homework1Cluster04.xml";
	private static final String TITLE = "Test matrices printed by GLIMMPSE for original Homework 1";
    private static final double[] ALPHA_LIST = {0.05};
    private static final double TOLERANCE = 0.001;

    private PowerChecker checker;

    @Before
	public void setUp() {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, true);
	}
	
    /**
     * Compare the calculated HLT approximate unconditional powers against simulation
     */
    @Test
    public void testPower()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();

        // add HLT as the statistical test
        params.addTest(GLMMTestFactory.Test.HOTELLING_LAWLEY_TRACE);

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int Q = 2;
        // create design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
        // add sample size multipliers
        //  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        params.addSampleSize(20);
        // build sigma matrix
        double CL1=15;
        double ICC1=0.13;
        double F1=(1 +(CL1-1)*ICC1)/CL1;
        double [][] sigma = {{1.21*F1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        params.addSigmaScale(1);

        // build beta matrix
        double [][] beta = {{1.24},{0.73}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        params.addBetaScale(1);

        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1, -1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        // build within subject contrast
        double [][] within = {{1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));


		checker.checkPower(params);
        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, true);
		assertTrue(checker.isSASDeviationBelowTolerance(TOLERANCE));
    }

}
