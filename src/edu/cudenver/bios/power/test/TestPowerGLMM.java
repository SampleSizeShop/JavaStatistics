package edu.cudenver.bios.power.test;

import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import edu.cudenver.bios.matrix.DesignEssenceMatrix;
import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.matrix.RandomColumnMetaData;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.Test;

import jsc.distributions.Normal;
import junit.framework.TestCase;

/**
 * Unit tests for power calculations on the GLMM.  Compares the calculated
 * values against a 10K iteration simulation using a t-test
 *
 */
public class TestPowerGLMM extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0.4997025,
    	0.8075886,
    	1.097641,
    	0.1651525,
    	0.2623301, 
    	0.3508015, 
    	0.1141548,
    	0.1812892,
    	0.2423835};
    private static final double[] SIGMA_SCALE_LIST = {1};
    private static final int[] SAMPLE_SIZE_LIST = {25};
    private Normal normalDist = new Normal();
    private DecimalFormat Number = new DecimalFormat("#0.000");

    /**
     * Test valid inputs for a univariate linear model with only fixed predictors
     */
    private void testValidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        System.out.println("Testing Univariate, Fixed");
        checkPower(params, true, true);
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
            List<Power> results = calc.getPower(params);
            fail();
        }
        catch (Exception e)
        {
            assertTrue(true);
        }
    }

    /**
     * Test case for a multivariate GLM with only fixed predictors
     */
    private void testValidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();

        System.out.println("Testing Multivariate, Fixed");
        checkPower(params, false, true);
    }


    /**
     * Test case for multivariate GLM with a baseline covariate.  Note that only
     * the Hotelling-Lawley and Univariate Repeated Measures tests are 
     * supported (see Glueck, Muller 2003 for details)
     */
    public void testValidMultivariateRandom()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        
        System.out.println("Testing Multivariate, Random, Quantile");
        checkPower(params, false, false);

    }


    /********** helper functions to create the matrices ***********/
    
    /**
     * Builds matrices for a univariate GLM with fixed predictors
     */
    private GLMMPowerParameters buildValidUnivariateInputs()
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
    private GLMMPowerParameters buildValidMultivariateFixedInputs()
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
    
    /**
     * Builds matrices for a multivariate GLM with a baseline covariate
     * Note, this matrix set matches the values produced in Table II from Glueck&Muller
     */
    private GLMMPowerParameters buildValidMultivariateRandomInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add power methods
        //for(PowerMethod method: PowerMethod.values()) params.addPowerMethod(method);
        //params.addPowerMethod(PowerMethod.CONDITIONAL_POWER);
        params.addPowerMethod(PowerMethod.QUANTILE_POWER);
        params.addQuantile(0.5);
        
        // add tests - only HL andUNIREP value for random case
        params.addTest(GLMMPowerParameters.Test.HOTELLING_LAWLEY_TRACE);
       // params.addTest(GLMMPowerParameters.Test.UNIREP);
        
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
        double [][] sigmaY = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

        // build sigma YG
        double [][] sigmaYG = {{0.5},{0.5}, {0.5}, {0}};
        params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));

        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0,0,0},{0,2,0,0},{0,0,0,0}};
        double [][] betaRandom = {{1,1,1,1}};
        params.setBeta(new FixedRandomMatrix(beta, betaRandom, false));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0,0,0},{0,0,0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0}, {1,0,-1}};
        double[][] betweenRandom = {{1}, {1}};
        params.setBetweenSubjectContrast(new FixedRandomMatrix(between, betweenRandom, true));
        
        // build within subject contrast
        double [][] within = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;     
    }
    
    /**
     * Run the power calculations for the specified parameters and tests
     * and assert whether they match simulation
     * 
     * @param params
     * @param testList
     */
    private void checkPower(GLMMPowerParameters params, 
            boolean univariate, boolean fixed)
    {
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();

        // build prefix string
        String uniStr = (univariate ? "U" : "M");
        String fixedStr = (fixed ? "F" : "R");

        int matches = 0;

        System.out.println("Calculating power...");
        List<Power> results = calc.getPower(params);

        System.out.println("Simulating power...");
        List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
        System.out.println("Multi?\tFixed?\tAlpha\tSigmaScale\tBetaScale\tTotal N\tPower\tPowerMethod\tQuantile");
        for(int i = 0; i < results.size() && i < simResults.size(); i++)
        {
            GLMMPower p = (GLMMPower) results.get(i);
            GLMMPower sp = (GLMMPower) simResults.get(i);
            System.out.println("Calculated: "+uniStr + "\t" + fixedStr + "\t" + 
                    p.getTest().toString() + "\t" + Number.format(p.getAlpha()) + "\t" + 
                    Number.format(p.getSigmaScale()) + "\t" + 
                    Number.format(p.getBetaScale()) + "\t" + 
                    p.getTotalSampleSize() + "\t" + 
                    Number.format(p.getActualPower()) + "\t" +
                    p.getPowerMethod() + "\t" + p.getQuantile());
            System.out.println("Simulated: "+uniStr + "\t" + fixedStr + "\t" + 
                    sp.getTest().toString() + "\t" + Number.format(sp.getAlpha()) + "\t" + 
                    Number.format(sp.getSigmaScale()) + "\t" + 
                    Number.format(sp.getBetaScale()) + "\t" + 
                    sp.getTotalSampleSize() + "\t" + 
                    Number.format(sp.getActualPower()) + "\t" +
                    p.getPowerMethod() + "\t" + p.getQuantile());
            if (powersAreSame(p.getActualPower(), sp.getActualPower())) matches++;
        }


        assertEquals(results.size(), matches);
    }
    
    /**
     * compares a calculated and simulated power by t-test and 
     * returns true if there is NO significant difference
     * 
     * @param calc calculated power value
     * @param sim simulated power value
     * @return true if the power do not differ significantly
     */
    private boolean powersAreSame(double calc, double sim)
    {
        double z = Math.abs((sim - calc) / Math.sqrt((sim * (1 - sim)) / SIMULATION_SIZE));
        
        double p = 2 * normalDist.upperTailProb(z);
        System.out.println("(p = " + Number.format(p) + ")");
        return (p > UNIT_TEST_ALPHA);
    }
}
