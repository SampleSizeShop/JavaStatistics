package edu.cudenver.bios.power.test.general;

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

public class TestSampleSizeGLMM extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] ALPHA_LIST = {0.05};    
    private static final double[] BETA_SCALE_LIST = {0,1,2};
    private static final double[] SIGMA_SCALE_LIST = {1,2};
    private static final double[] POWER_LIST = {0.7,0.8,0.9};
    private Normal normalDist;
    private DecimalFormat Number;

    public void setUp()
    {
        normalDist = new Normal();
        Number = new DecimalFormat("#0.000");

    }
    
    private void testValidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
                
        List<Power> results = calc.getSampleSize(params);
        //List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
        System.out.println("Multi?\tFixed?\tAlpha\tSigmaScale\tBetaScale\tTotal N\tPower");
        for(Power power: results)
        {
        	GLMMPower p = (GLMMPower) power;
        	System.out.println("U\tF\t" + p.getTest().toString() + "\t" + Number.format(p.getAlpha()) +
        			"\t" + Number.format(p.getSigmaScale()) + "\t" + Number.format(p.getBetaScale()) + 
        			"\t" + p.getTotalSampleSize() + "\t" + Number.format(p.getActualPower()));
        }
    }

    private void testInvalidUnivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidUnivariateInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);

        try
        {
        	List<Power> results = calc.getSampleSize(params);
        	//List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
        	for(Power p: results)
        	{
        		System.out.println("Univariate, Fixed: " + p.toXML());
        	}
        }
        catch (Exception e)
        {
            System.out.println("Exception: " + e.getMessage());
        }
    }

    public void testValidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        
        List<Power> results = calc.getSampleSize(params);
        System.out.println("Multi?\tFixed?\tAlpha\tSigmaScale\tBetaScale\tTotal N\tPower");
        for(Power power: results)
        {
        	GLMMPower p = (GLMMPower) power;
        	System.out.println("M\tF\t" + p.getTest().toString() + "\t" + Number.format(p.getAlpha()) +
        			"\t" + Number.format(p.getSigmaScale()) + "\t" + Number.format(p.getBetaScale()) + 
        			"\t" + p.getTotalSampleSize() + "\t" + Number.format(p.getActualPower()));
        }

    }

    private void testInvalidMultivariateFixed()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateFixedInputs();
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        params.setBeta(null);

        try
        {
                List<Power> results = calc.getSampleSize(params);
                //List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
                for(Power p: results)
                {
                    System.out.println("Multivariate, Fixed: " + p.toXML());
                }
        }
        catch (Exception e)
        {
            System.out.println("Exception: " + e.getMessage());
        }
    }

    private void testValidMultivariateRandom()
    {
        // build the inputs
        GLMMPowerParameters params = buildValidMultivariateRandomInputs();
        
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        
        List<Power> results = calc.getSampleSize(params);
        //List<Power> simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
        for(Power p: results)
        {
        	System.out.println("Multivariate, Random: " + p.toXML());
        }
    }

    private void testInvalidMultivariateRandom()
    {

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
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
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
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
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
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
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
    
    private boolean powersAreSame(double calc, double sim)
    {
        double z = Math.abs((sim - calc) / Math.sqrt((sim * (1 - sim)) / SIMULATION_SIZE));
        
        double p = 2 * normalDist.upperTailProb(z);
        System.out.println("(p = " + Number.format(p) + ")");
        return (p > UNIT_TEST_ALPHA);
    } 
    
}

