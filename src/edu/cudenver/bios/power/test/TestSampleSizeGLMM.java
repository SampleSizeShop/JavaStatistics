package edu.cudenver.bios.power.test;

import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

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
        params.setPowerMethod(PowerMethod.QUANTILE_POWER);
        params.setQuantile(0.25);
        
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


    private GLMMPowerParameters buildValidUnivariateInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
       
        // add tests
        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
        {
            if (test != GLMMPowerParameters.Test.NONE) params.addTest(test);
        }
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new Array2DRowRealMatrix(beta));
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
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(10,1));
        params.setDesignEssence(essenceMatrix);
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));

        return params;     
    }   


    private GLMMPowerParameters buildValidMultivariateFixedInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests
        for(GLMMPowerParameters.Test test: GLMMPowerParameters.Test.values()) 
        {
            if (test != GLMMPowerParameters.Test.NONE) params.addTest(test);
        }
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int Q = 4;
        // create design matrix
        RealMatrix essenceData = MatrixUtils.createRealIdentityMatrix(Q);
        EssenceMatrix essence = new EssenceMatrix(essenceData);
        essence.setRowMetaData(0, new RowMetaData(5,1));
        essence.setRowMetaData(1, new RowMetaData(5,1));
        essence.setRowMetaData(2, new RowMetaData(5,1));
        essence.setRowMetaData(3, new RowMetaData(5,1));
        params.setDesignEssence(essence);
        // add powers
        for(double power: POWER_LIST) params.addPower(power);
        
        // build sigma matrix
        double rho = 0.4;
        double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}}; // compound symmetry
        //double [][] sigma = {{1,0.2,0.3},{0.2,1,0.2},{0.3,0.2,1}}; // toeplitz
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // add sigma scale values
        for(double sigmaScale: SIGMA_SCALE_LIST) params.addSigmaScale(sigmaScale);
        
        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));

        // build within subject contrast
        double [][] within = {{1,1},{-1,0},{0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        //RealMatrix U = params.getWithinSubjectContrast();
        //RealMatrix upu = U.multiply(U.transpose());
        
        return params;     
    }   
    
    private boolean powersAreSame(double calc, double sim)
    {
        double z = Math.abs((sim - calc) / Math.sqrt((sim * (1 - sim)) / SIMULATION_SIZE));
        
        double p = 2 * normalDist.upperTailProb(z);
        System.out.println("(p = " + Number.format(p) + ")");
        return (p > UNIT_TEST_ALPHA);
    }
    
    private GLMMPowerParameters buildValidMultivariateRandomInputs()
    {
        GLMMPowerParameters params = new GLMMPowerParameters();
        
        // add tests - only HL andUNIREP value for random case
        params.addTest(GLMMPowerParameters.Test.HOTELLING_LAWLEY_TRACE);
        params.addTest(GLMMPowerParameters.Test.UNIREP);
        
        // add alpha values
        for(double alpha: ALPHA_LIST) params.addAlpha(alpha);

        int P = 3;
        int Q = 3;
        // create design matrix
        double[][] essData = {{1,0,0},{0,1,0}};
        RealMatrix essenceData = new Array2DRowRealMatrix(essData);
        EssenceMatrix essence = new EssenceMatrix(essenceData);
        essence.setRowMetaData(0, new RowMetaData(5,1));
        essence.setRowMetaData(1, new RowMetaData(5,1));
        //essence.setRowMetaData(2, new RowMetaData(5,1));
        //essence.setRowMetaData(3, new RowMetaData(5,1));
        params.setDesignEssence(essence);
        
        // add a random column with specified mean/variance
        ColumnMetaData randColMD = new ColumnMetaData();
        randColMD.setPredictorType(PredictorType.RANDOM);
        randColMD.setMean(MEAN);
        randColMD.setVariance(VARIANCE);
        essence.setColumnMetaData(2, randColMD);
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
        double [][] beta = {{0.25,0},{0,0.25},{0.9,0}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // add beta scale values
        for(double betaScale: BETA_SCALE_LIST) params.addBetaScale(betaScale);
        
        // build theta null matrix
        double [][] theta0 = {{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,0,0}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
        
        // build within subject contrast
        double [][] within = {{1,0},{0,1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;     
    }
    
    
}

