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

import java.io.File;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.test.PowerChecker;
import junit.framework.TestCase;

/**
 * Unit test for fixed univariate design with comparison against
 * simulation and SAS output.
 * @author Sarah Kreidler
 *
 */
public class TestConditionalMultivariateInteraction extends TestCase
{
	private static final String DATA_FILE =  "sas" + File.separator + "data" + File.separator + "TestConditionalMultivariateInteraction.xml";
	private static final String OUTPUT_FILE = "text" + File.separator + "results" + File.separator + "TestConditionalMultivariateInteraction.html";
	private static final String TITLE = "Power results for multivariate interaction";
	private PowerChecker checker;
	
	public void setUp()
	{
		try
		{
			checker = new PowerChecker(DATA_FILE, true);
		}
		catch (Exception e)
		{
			System.err.println("Setup failed: " + e.getMessage());
			fail();
		}
	}
	
    /**
     * Test valid inputs for a univariate linear model with only fixed predictors
     */
    public void testMultivariateInteraction()
    {
        // build the inputs
        GLMMPowerParameters params = new GLMMPowerParameters();
        	
        // build the matrix inputs
        
        // add tests
        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
        {
            params.addTest(test);
        }
        
        // add alpha values
        params.addAlpha(0.01);

        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        for(double scale = 0; scale <= 2.0; scale += 0.50) params.addBetaScale(scale);
        
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
        double[][] essenceData = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        RowMetaData[] rowMd = {new RowMetaData(1), new RowMetaData(1), 
        		new RowMetaData(1), new RowMetaData(1)};
        DesignEssenceMatrix essenceMatrix = new DesignEssenceMatrix(essenceData, rowMd, null, null);
        params.setDesignEssence(essenceMatrix);
        // add sample size multipliers
        params.addSampleSize(5);
        params.addSampleSize(10);
        
        // build between subject contrast
        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));
        
        // build within subject contrast
        double [][] within = {{1,1},{-1,0},{0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));
        	
        System.out.println(TITLE);
        checker.checkPower(params);
		checker.outputResults();
		checker.outputResults(TITLE, OUTPUT_FILE);
		assertTrue(checker.isSASDeviationBelowTolerance());
		assertTrue(checker.isSimulationDeviationBelowTolerance());
		checker.reset();
    }

    /**
     * Tests if the calculator throws an exception on invalid inputs
     */
    public void testInvalidInputs()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);
        
        try
        {
        	calc.getPower(params);
            fail();
        }
        catch (Exception e)
        {
            assertTrue(true);
        }
    }

    


    /********** helper functions to create the matrices ***********/
    
    /**
     * Builds matrices for a univariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildValidUnivariateInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
       


        return params;     
    }   
}
