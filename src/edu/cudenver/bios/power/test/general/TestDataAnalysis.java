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

import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.Test;
import jsc.distributions.FishersF;
import jsc.distributions.Normal;
import junit.framework.TestCase;

/**
 * Test the GLMMTest objects observed F, degrees of freedom when
 * doing data analysis - used for simulation.  This isn't really an
 * automated unit test - used to compare manually against LINMOD 
 * library (SAS/IML) written by Keith Muller.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestDataAnalysis extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {2};
    private static final double[] SIGMA_SCALE_LIST = {2};
    private static final int[] SAMPLE_SIZE_LIST = {5};
    private Normal normalDist = new Normal();
    private DecimalFormat number = new DecimalFormat("#0.000");
    
    private double[][] univariateErrorData = {
    	{ 2.67 },
    	{ 0.97 },
    	{0.55},
    	{1.57},
    	{-0.12},
    	{-1.02},
    	{0.94},
    	{-1.25},
    	{0.94},
    	{-0.08},
    	{0.09},
    	{0.78},
    	{-0.30},
    	{-0.59},
    	{-0.53},
    	{0.20},
    	{2.89},
    	{-1.04},
    	{2.05},
    	{-1.08}
    };
    private Array2DRowRealMatrix univariateErrorMatrix = new Array2DRowRealMatrix(univariateErrorData);
    
    private double[][] multivariateErrorData = {
    		{-1.56,-0.88, 0.58},
    		 {2.03, -1.88, -0.16},
    		 {0.19, -0.07,  0.67},
    		 {0.83, -0.22,  0.63},
    		 {-1.40, -1.47,  1.71},
    		 {0.07, -0.53, -0.09},
    		 {-0.20, -0.87,  0.06},
    		 {0.27, -1.12, -1.42},
    		 {-0.01,  0.08, -0.38},
    		{1.66, -0.53, -0.75},
    		{0.25,  1.73,  0.08},
    		{-1.05,  0.61, -0.15},
    		{2.01, -0.54, -1.34},
    		{-1.06,  0.16,  1.88},
    		{-1.43,  2.18, -0.39},
    		{-0.75, -0.19, -1.46},
    		{-1.08, -0.91,  0.80},
    		{-1.16, -1.68,  0.61},
    		{1.67,  1.00,  0.55},
    		{-0.47, -0.07,  0.14}
        };
        private Array2DRowRealMatrix multivariateErrorMatrix = new Array2DRowRealMatrix(multivariateErrorData);
    
    
    private class FitResult {
    	public double ndf;
    	public double ddf;
    	public double fCrit;
    	public double pValue;
    	
    	public FitResult(double ndf, double ddf, double fCrit, double pValue)
    	{
    		this.ndf = ndf;
    		this.ddf = ddf;
    		this.fCrit = fCrit;
    		this.pValue = pValue;
    	}
    };  
    
	private void testUnirep()
	{
		GLMMPowerParameters params = buildUnivariateInputs();
		params.addTest(Test.UNIREP);
		FitResult result = fitModel(params, univariateErrorMatrix);

		
	}
	
	private void testBox()
	{
		GLMMPowerParameters params = buildUnivariateInputs();
		params.addTest(Test.UNIREP_BOX);
		FitResult result = fitModel(params, univariateErrorMatrix);

	}
	
	private void testGeisserGreenhouse()
	{
		GLMMPowerParameters params = buildUnivariateInputs();
		params.addTest(Test.UNIREP_GEISSER_GREENHOUSE);
		FitResult result = fitModel(params, univariateErrorMatrix);
		
	}
	
	private void testHuynhFeldt()
	{
		GLMMPowerParameters params = buildUnivariateInputs();
		params.addTest(Test.UNIREP_HUYNH_FELDT);
		FitResult result = fitModel(params, univariateErrorMatrix);
		
	}
	
	public void testWilksLambda()
	{
//		GLMMPowerParameters paramsUni = buildUnivariateInputs();
//		
//		paramsUni.addTest(Test.WILKS_LAMBDA);
//		FitResult resultUni = fitModel(paramsUni, univariateErrorMatrix);
		
		GLMMPowerParameters params = buildMultivariateInputs();
		
		params.addTest(Test.WILKS_LAMBDA);
		FitResult result = fitModel(params, multivariateErrorMatrix);

	}
	
	public void testPillaiBartlett()
	{
		
		GLMMPowerParameters params = buildMultivariateInputs();
		
		params.addTest(Test.PILLAI_BARTLETT_TRACE);
		FitResult result = fitModel(params, multivariateErrorMatrix);

	}
	
	public void testHotellingLawley()
	{
		GLMMPowerParameters params = buildMultivariateInputs();
		
		params.addTest(Test.HOTELLING_LAWLEY_TRACE);
		FitResult result = fitModel(params, multivariateErrorMatrix);

	}
	
    /**
     * Builds matrices for a univariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildUnivariateInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
       
        // add tests
        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
        {
            params.addTest(test);
        }
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new FixedRandomMatrix(beta, null, false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build design matrix
        double[][] essenceData = {{1,0},{0,1}};
        RowMetaData[] rowMd = {new RowMetaData(1), new RowMetaData(1)};
        DesignEssenceMatrix essenceMatrix = new DesignEssenceMatrix(essenceData, rowMd, null, null);
        params.setDesignEssence(essenceMatrix);
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        return params;     
    }   

    /**
     * Builds matrices for a multivariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildMultivariateInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
     
        // add tests
//        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
//        {
//            params.addTest(test);
//        }
        params.addTest(Test.WILKS_LAMBDA);
        params.addTest(Test.PILLAI_BARTLETT_TRACE);
        params.addTest(Test.HOTELLING_LAWLEY_TRACE);

        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int Q = 4;
        // create design matrix
        RealMatrix essenceData = MatrixUtils.createRealIdentityMatrix(Q);
        RowMetaData[] rowMd = {
        		new RowMetaData(1), 
        		new RowMetaData(1), 
        		new RowMetaData(1), 
        		new RowMetaData(1)
        		};
        DesignEssenceMatrix essence = 
        	new DesignEssenceMatrix(essenceData.getData(), rowMd, null, null);
        params.setDesignEssence(essence);
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build sigma matrix
        double rho = 0.4;
        double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}}; // compound symmetry
        // double [][] sigma = {{1,0.2,0.3},{0.2,1,0.2},{0.3,0.2,1}}; // toeplitz
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new FixedRandomMatrix(beta, null,  false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, null, true));

        // build within subject contrast
        double [][] within = {{1,1},{-1,0},{0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        //RealMatrix U = params.getWithinSubjectContrast();
        //RealMatrix upu = U.multiply(U.transpose());
        
        return params;     
    }   
	
    private FitResult fitModel(GLMMPowerParameters params, RealMatrix error)
    {
    	params.getFirstAlpha();
    	params.getFirstBetaScale();
    	params.getFirstSigmaScale();
    	params.getFirstSampleSize();
    	params.getFirstTest();
    	
        RealMatrix scaledBeta = params.getScaledBeta();
        RealMatrix scaledSigma = params.getScaledSigmaError();
        RealMatrix X = params.getDesign();
        int N = X.getRowDimension();
        int rankX = new SingularValueDecompositionImpl(X).getRank();
        
        // calculate simulated Y based on Y = X beta + error
        //RealMatrix Ysim = (X.multiply(scaledBeta)).add(error);
        double[][] simData = {
            	{1,2,3},
            	{1,2,3},
            	{1,2.7, 3},
            	{1,1,3},
            	{1,2, 3},
            	{1,4, 2},
            	{1,2, 2},
            	{1,3.5, 2},
            	{1,4, 2},
            	{1,4, 2},
            	{1,4, 2},
            	{1,6 ,2}, 
            	{1,4, 2}, 
            	{1,5 ,2}, 
            	{1,4, 2}, 
            	{-1,0, 2}, 
            	{-1,0, 1}, 
            	{-1,0, 2}, 
            	{-1.2, 0, 1.5}, 
            	{-1,0, 2}
            };
        RealMatrix Ysim = new Array2DRowRealMatrix(simData);
        
//        for(int r = 0; r < Ysim.getRowDimension(); r++)
//        {
//        	for (int c = 0; c < Ysim.getColumnDimension(); c++)
//        	{
//        		System.out.print(Ysim.getEntry(r, c) + " ");
//        	}
//        	System.out.print("\n");
//        }
        // calculate beta-Hat
        RealMatrix XtX = X.transpose().multiply(X);
        RealMatrix XtXInverse = new LUDecompositionImpl(XtX).getSolver().getInverse();
        RealMatrix betaHat = XtXInverse.multiply(X.transpose()).multiply(Ysim);
//      for(int r = 0; r < betaHat.getRowDimension(); r++)
//      {
//      	for (int c = 0; c < betaHat.getColumnDimension(); c++)
//      	{
//      		System.out.print(betaHat.getEntry(r, c) + " ");
//      	}
//      	System.out.print("\n");
//      }
        // calculate Y-hat
        RealMatrix YHat = (X.multiply(betaHat));
        // calculate sigma-Hat
        RealMatrix Ydiff = Ysim.subtract(YHat);
        RealMatrix sigmaHat = (Ydiff.transpose().multiply(Ydiff)).scalarMultiply(((double) 1/(double)(N - rankX)));    
        
        // set the scaled sigma / beta to the new values
        params.setScaledBeta(betaHat);
        params.setScaledSigmaError(sigmaHat);
        
        // calculate the observed F for the simulation
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        double fobs = glmmTest.getObservedF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

        // get the p-value from a central F distribution
        double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
        double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
        
        FishersF fdist = new FishersF(ndf, ddf);
        double pvalue = 1 - fdist.cdf(fobs);
        
		System.out.println("Test: " + params.getCurrentTest() + 
				" Ndf: " + number.format(ndf) + 
				" Ddf: " + number.format(ddf) +
				" F-crit: " + number.format(fobs) + 
				" p-value: " + number.format(pvalue));
		
        return new FitResult(ndf, ddf, fobs, pvalue);
    }
}
