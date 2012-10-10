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

import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import junit.framework.TestCase;

public class TestLawANOVA extends TestCase
{
    private static DecimalFormat Number = new DecimalFormat("#0.0000");

	/**
	 * Reproduce the list of power displayed in the figure associated
	 * with the Law ANCOVA example in the GLIMMPSE manuscript
	 */
	public void testANOVA()
	{
		try
		{
	        // build the inputs
	        GLMMPowerParameters params = new GLMMPowerParameters();
	        	
	        // add alpha values
	        params.addAlpha(0.01);
	        
	        // build the design matrix - reference cell coded
	        params.setDesignEssence(MatrixUtils.createRealIdentityMatrix(4));
	        
	        // add desired power
	        params.addPower(0.9);
	        
	        // build beta matrix
	        double [][] beta = {{5},{0},{0},{0}};
	        params.setBeta(new FixedRandomMatrix(beta, null, false));
	        // add beta scale values
	        params.addBetaScale(1);
	        
	        // build between subject contrast
	        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
	        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));
	        
	        // build within subject contrast
	        double [][] within = {{1}};
	        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));
	        	
	        // build theta null matrix
	        double [][] theta0 = {{0},{0},{0}};
	        params.setTheta(new Array2DRowRealMatrix(theta0));
	        
	        // build covariance of E (SD in manuscript = 8.1 => variance = 65.61
	        double [][] sigmaE = {{65.61}}; 
	        params.setSigmaError(new Array2DRowRealMatrix(sigmaE));
	        // add sigma scale factors
	        params.addSigmaScale(1);
	        params.addSigmaScale(0.5);
	        params.addSigmaScale(2);
	        
	        // add tests
	        params.addTest(Test.UNIREP);

	        // get the power results
	        GLMMPowerCalculator calc = new GLMMPowerCalculator();
	        List<Power> results = calc.getSampleSize(params);
	        
	        // output the sample size results
			writeResults(results);
			assertTrue(true);			
		}
		catch (Exception e)
		{
			System.err.print(e.getMessage());
			fail();
		}
	}
	
	/**
	 * 
	 * @param results list of power results
	 */
	private void writeResults(List<Power> results)
	{
		System.out.println("Test, Actual Power, Total Sample Size," +
				"Beta Scale, Sigma Scale, Alpha, Nominal Power, Power Method");

		for(Power power: results)
		{
			GLMMPower glmmPower = (GLMMPower) power;
			System.out.println(testToString(glmmPower.getTest()) + ", " + 
					Number.format(glmmPower.getActualPower()) + ", " + 
					glmmPower.getTotalSampleSize() + ", " + 
					glmmPower.getBetaScale() + ", " + 
					glmmPower.getSigmaScale() + ", " + 
					Number.format(glmmPower.getNominalPower()) + ", " + 
					powerMethodToString(glmmPower.getPowerMethod()));
		}
	}
	
	/**
	 * Pretty print test name
	 */
	private String testToString(Test test)
	{
		switch (test)
		{
			case HOTELLING_LAWLEY_TRACE:
				return "hlt";
			case WILKS_LAMBDA:
				return "wl";
			case PILLAI_BARTLETT_TRACE:
				return "pbt";
			case UNIREP:
				return "unirep";
			case UNIREP_GEISSER_GREENHOUSE:
				return "unirepGG";
			case UNIREP_BOX:
				return "unirepBox";
			case UNIREP_HUYNH_FELDT:
				return "unirepHF";				
		}
		return "";
	}
	
	/**
	 * Pretty print power method
	 */
	private String powerMethodToString(PowerMethod method)
	{
		switch (method)
		{
		case CONDITIONAL_POWER:
			return "conditional";
		case UNCONDITIONAL_POWER:
			return "unconditional";
		case QUANTILE_POWER:
			return "quantile";
		}
		return "";
	}
}
