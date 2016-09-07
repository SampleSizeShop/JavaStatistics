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
import edu.cudenver.bios.matrix.OrthogonalPolynomials;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import edu.cudenver.bios.power.test.ValidationReportBuilder;
import edu.cudenver.bios.utils.Factor;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * Unit test for polynomial trends in U matrix
 * with comparison against simulation and SAS output.
 * Based on example 7 from POWERLIB (Johnson et al., 2009)
 * @author Sarah Kreidler
 *
 */
public class ConditionalOrthogonalPolynomial1FactorTest {

    private static final String DATA_FILE =  "TestConditionalOrthogonalPolynomial1Factor.xml";
    private static final String TITLE = "GLMM(F) Example 7. Power for a time by " +
            "treatment interaction using orthogonal polynomial contrast for time";
    private static final double TOLERANCE = 0.000001;

    private PowerChecker checker;
    private boolean verbose = false;

    @Before
    public void setUp() {
        List<GLMMPower> sasPowers = Utils.readSasPowers(DATA_FILE);
        checker = new PowerChecker(sasPowers, false);
    }

    /**
     * Test GLMM(F) with polynomial contrasts in U matrix
     */
    @Test
    public void testPolynomial1Factor() {
        // build the inputs
        GLMMPowerParameters params = new GLMMPowerParameters();

        // add tests
        for(GLMMTestFactory.Test test: GLMMTestFactory.Test.values()) {
            params.addTest(test);
        }

        // add alpha values - bonferroni corrected for 6 comparisons
        params.addAlpha(0.05);

        // build beta matrix
        double [][] beta = {{0,0,0,0,1},{1,0,0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        params.addBetaScale(1);

        // build theta null matrix
        double [][] theta0 = {{0,0,0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build sigma matrix
        double rho = 0.375;
        double var = 1.5;
        double [][] sigma = {
                {var,rho,rho,rho,rho},
                {rho,var,rho,rho,rho},
                {rho,rho,var,rho,rho},
                {rho,rho,rho,var,rho},
                {rho,rho,rho,rho,var}
        };
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        params.addSigmaScale(1);

        // build design matrix
        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(2));
        // add sample size multipliers
        params.addSampleSize(10);
        params.addSampleSize(20);
        params.addSampleSize(40);

        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        // build within subject contrast
        double[] times ={2, 4 ,6, 8, 10};
        String name = "times";
        ArrayList<Factor> factorList = new ArrayList<Factor>();
        Factor timeFactor = new Factor(name, times);
        factorList.add(timeFactor);
        RealMatrix U =
            OrthogonalPolynomials.withinSubjectContrast(factorList).getMainEffectContrast(timeFactor).getContrastMatrix();
        if (verbose) printMatrix("U Matrix", U);
        params.setWithinSubjectContrast(U);

        checker.checkPower(params);

        ValidationReportBuilder reportBuilder = new ValidationReportBuilder();
        reportBuilder.createValidationReportAsStdout(checker, TITLE, false);
        assertTrue("SAS deviation " + checker.getMaxSasDeviation() + " is not below tolerance " + TOLERANCE,
                    checker.isSASDeviationBelowTolerance(TOLERANCE));
    }

    /**
     * Write the matrix to std out
     * @param m
     */
    private void printMatrix(String title, RealMatrix m)
    {
        System.out.println(title);
        DecimalFormat Number = new DecimalFormat("#0.000");
        for(int row = 0; row < m.getRowDimension(); row++)
        {
            for(int col= 0; col < m.getColumnDimension(); col++)
            {
                System.out.print(Number.format(m.getEntry(row, col)) + "\t");
            }
            System.out.print("\n");
        }
    }
}
