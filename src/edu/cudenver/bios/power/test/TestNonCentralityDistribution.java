package edu.cudenver.bios.power.test;

import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;

import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RandomColumnMetaData;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.glmm.NonCentralityDistribution;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import jsc.distributions.Normal;
import junit.framework.TestCase;

public class TestNonCentralityDistribution extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0,0.5,1,1.5,2};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final int[] SAMPLE_SIZE_LIST = {5};
    private Normal normalDist = new Normal();
    private DecimalFormat Number = new DecimalFormat("#0.000");
    
    
    public void testApproximateNonCentralCDF()
    {
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        params.getFirstAlpha();
        params.getFirstBetaScale();
        params.getFirstSigmaScale();
        params.getFirstSampleSize();
        params.getFirstTest();
        
        NonCentralityDistribution ncd = 
            new NonCentralityDistribution(params, false);

        for(double w = 0.10; w < 1; w += 0.1)
        {
            double nonCentralityParam = ncd.inverseCDF(w);
            System.out.println("Quantile: " + w + " inverseCDF: " + nonCentralityParam);
        }
        
    }
    
    /**
     * Builds matrices for a multivariate GLM with a baseline covariate
     */
    private GLMMPowerParameters buildValidMultivariateRandomInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add power methods
        //for(PowerMethod method: PowerMethod.values()) params.addPowerMethod(method);
        params.addPowerMethod(PowerMethod.CONDITIONAL_POWER);
        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
        params.addQuantile(0.25);
        params.addQuantile(0.5);
        params.addQuantile(0.75);
        
        // add tests - only HL andUNIREP value for random case
        params.addTest(GLMMPowerParameters.Test.HOTELLING_LAWLEY_TRACE);
        params.addTest(GLMMPowerParameters.Test.UNIREP);
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int P = 3;
        int Q = 3;
        // create design matrix
        double[][] essFixedData = {{1,0,0},{0,1,0},{0,0,1}};
        RowMetaData[] rowMd = {
        		new RowMetaData(1), 
        		new RowMetaData(1), 
        		new RowMetaData(1)
        		};
        double[][] essRandomData = {{1},{1},{1}};
        RandomColumnMetaData[] randColMd = {new RandomColumnMetaData(MEAN, VARIANCE)};
        DesignEssenceMatrix essence = new DesignEssenceMatrix(essFixedData, rowMd, 
        		essRandomData, randColMd);
        params.setDesignEssence(essence);
        // add sample size multipliers
        for(int sampleSize: SAMPLE_SIZE_LIST) params.addSampleSize(sampleSize);
        
        // build sigma G matrix
        double[][] sigmaG = {{1}};
        params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));

        // build sigma Y matrix
        double rho = 0.4;
        double [][] sigmaY = {{1,0},{0,1}};
        params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

        // build sigma YG
        double rhoYG = 0.8;
        double [][] sigmaYG = {{0.9},{0}};
        params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0},{0,0},{0,0}};
        double [][] betaRandom = {{1},{1}};
        params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,0,0}};
        double[][] betweenRandom = {{1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));
        
        // build within subject contrast
        double [][] within = {{1,0},{0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;     
    }
    
}
