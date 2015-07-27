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

import java.io.File;
import java.util.List;

import edu.cudenver.bios.power.GLMMPower;
import junit.framework.TestCase;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Unit test for fixed multvariate design with comparison against
 * simulation and SAS output.  Based on example 5 from powerlib 
 * (Johnson et al., 2009)
 * @author Sarah Kreidler
 *
 */
public class ConditionalMultivariateInteractionTest {

	private static final String DATA_FILE =  "TestConditionalMultivariateInteraction.xml";
	private static final String TITLE =
	        "GLMM(F) Example 5. Power for a test of interaction in a multivariate model";
    private static final double TOLERANCE = 0.1;

	private PowerChecker checker;

    @Before
	public void setUp() {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, false);
	}
	
    /**
     * Test valid inputs for a multivariate GLMM(F) testing for interaction effect
     */
    @Test
    public void testMultivariateInteraction()
    {
        // build the inputs
        GLMMPowerParameters params = new GLMMPowerParameters();
        	
        // build the matrix inputs
        
        // add tests
        for(GLMMTestFactory.Test test: GLMMTestFactory.Test.values()) {
            params.addTest(test);
        }

        // add alpha values
        params.addAlpha(0.01);

        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        for(double scale = 0; scale <= 2.0; scale += 0.50) {
            params.addBetaScale(scale);
        }
        
        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double rho = 0.4;
        double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        params.addSigmaScale(1);
        params.addSigmaScale(2);
        
        // build design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(4));
        // add sample size multipliers
        params.addSampleSize(5);
        params.addSampleSize(10);
        
        // build between subject contrast
        double [][] between = {
                {1,-1, 0, 0},
                {1, 0,-1, 0},
                {1, 0, 0,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));
        
        // build within subject contrast
        double [][] within = {
                { 1, 1},
                {-1, 0},
                { 0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        checker.checkPower(params);

        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, false);

        assertTrue("results outside tolerance: " + TOLERANCE, checker.isSASDeviationBelowTolerance(TOLERANCE));
    }
}
