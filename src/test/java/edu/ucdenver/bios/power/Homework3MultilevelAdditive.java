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
public class Homework3MultilevelAdditive {

	private static final String DATA_FILE =  "Homework3MultilevelAdditive.xml";
	private static final String TITLE = "Homework03Additive.SAS Test matrices printed by GLIMMPSE for Homework";
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

        int Q = 3;
        // create design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(Q));
        // add sample size multipliers
        //  for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        params.addSampleSize(15);
        // build sigma matrix
        double [][] sigma = {{19.360, 16.632, 0.528} ,
                             {16.632, 17.640, 1.008} ,
                             {0.528,  1.008, 0.360 }};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));

        double [][] Rho = {{0.04},{0.07}};
        double [][] K = {{5},{4}};
        params.addSigmaScale(calculateSigScale(Rho, K));

        // build beta matrix
        double [][] beta = {{.3,.3, .3},
                            {.1,.1, .1},
                            {.1,.1, .1 }};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        params.addBetaScale(1);

        // build theta null matrix
        double [][] theta0 = {{0, 0, 0},
                              {0, 0, 0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{-1, 1, 0},
                               {-1, 0, 1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        // build within subject contrast
        params.setWithinSubjectContrast(MatrixUtils.createRealIdentityMatrix(3));


		checker.checkPower(params);
        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, true);
		assertTrue(checker.isSASDeviationBelowTolerance(TOLERANCE));
    }

    private double calculateSigScale(double [][] Rho, double [][]K) {
        //*Inputs: RHO= D x 1 vector of correlations in interval [0,1); <- column
        //*        K  = D x 1 vector of nested cluster dimensions; <-
        //*Assumes integer D>0;  *Limit D to <=10 or so;
        //*Level 1 nested within level 2, nested within level 3, etc.;
        double FACTOR=1;
        double KSTARD=1;
        int DMAX=Rho.length;
        for (int D=0; D < DMAX; D++) {
            KSTARD = KSTARD * K[D][0];
            FACTOR = FACTOR / K[D][0] + Rho[D][0] * (1 - 1 / KSTARD);
        }
        return FACTOR;
    }

}
