package edu.cudenver.bios.powersamplesize.test;

import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.PowerGLMM;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.PowerMethod;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;

import junit.framework.TestCase;

public class TestPowerGLMM extends TestCase
{
    private static final int SIMULATION_SIZE = 100000;
    private static final double PRECISION = 0.01;
    private static final double ALPHA = 0.01;    
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 4.3;
    private static final double[] BETA_SCALE = {0, 0.5,2};
    private static final double[] SIGMA_SCALE = {1,2};
    private static final int[] SAMPLE_SIZE = {20,40};

    private void testValidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();

        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPower("Valid Univariate, Fixed, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Univariate, Fixed, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Univariate, Fixed, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Univariate, Fixed, W",calc, goodParams);        
    }

    private void testInvalidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();
        goodParams.setBeta(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Univariate, Fixed, UNIREP", calc, goodParams);

    }

    private void testValidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs(false);

        PowerGLMM calc = new PowerGLMM();
//        goodParams.setTestStatistic(TestStatistic.UNIREP);
//        checkPower("Valid Multivariate, Fixed, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Multivariate, Fixed, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Multivariate, Fixed, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Multivariate, Fixed, W",calc, goodParams);        

    }

    private void testInvalidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs(false);

        goodParams.setSigmaError(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Multivariate, Fixed, UNIREP", calc, goodParams);

    }

    public void testValidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs(true);
        goodParams.setPowerMethod(PowerMethod.QUANTILE_POWER);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPower("Valid Multivariate, Random, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Multivariate, Random, HLT",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.PILLAU_BARTLETT_TRACE);
        checkPower("Valid Multivariate, Random, PB",calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
        checkPower("Valid Multivariate, Random, W",calc, goodParams);  



    }

    public void testInvalidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateInputs(true);



    }

    private void checkPowerFail(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        try
        {
            double calculated = calc.getCalculatedPower(params);
            double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) params, SIMULATION_SIZE);
            assert(!(Math.abs(simulated - calculated) < PRECISION));
        }
        catch(Exception e)
        {
            System.out.println(label + ">> Recognized invalid inputs");
            assert(true);
        }

    }

    private void checkPower(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        int tests = 0; // number of tests run
        int matches = 0; // number of matches between calculated and simulated power
        DecimalFormat Number = new DecimalFormat("#0.0000");


        for(double sigmaScale : SIGMA_SCALE)
        {
            for(double betaScale : BETA_SCALE)
            {
                for(int sampleSize: SAMPLE_SIZE)
                {
                    tests++;     
                    LinearModelPowerSampleSizeParameters testParams = 
                        new LinearModelPowerSampleSizeParameters(params);
                    // scale the inputs
                    testParams.setBeta(params.getBeta().scalarMultiply(betaScale));
                    //testParams.setSigmaError(params.getSigmaError().scalarMultiply(sigmaScale));
                    testParams.setSampleSize(sampleSize);
                    try
                    {
                        double calculated = calc.getCalculatedPower(testParams);
                        double simulated = calc.getSimulatedPower(testParams, SIMULATION_SIZE);

                        System.out.println(label + "["+sigmaScale+","+betaScale+","+sampleSize
                                +"]>> Calculated power: " + Number.format(calculated) + ", simulated power: " + Number.format(simulated));
                        if (Math.abs(simulated - calculated) < PRECISION) matches++;
                    }
                    catch(Exception e)
                    {
                        System.out.println(label + ">> Failed to calculate power: " + e.getMessage());
                    }        

                }
            }
        }
        assertEquals(tests, matches);
    }

    private LinearModelPowerSampleSizeParameters buildValidUnivariateInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);

        // build beta matrix
        double [][] beta = {{0},{1}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // build design matrix
        double[][] essenceData = {{1,0},{0,1}};
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(10,1));
        params.setDesignEssence(essenceMatrix);

        // build between subject contrast
        double [][] between = {{1,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));

        return params;     
    }   


    private LinearModelPowerSampleSizeParameters buildValidMultivariateInputs(boolean hasRandom)
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);

        int P = 3;
        int Q = 4;
        // create design matrix
        RealMatrix essenceData = MatrixUtils.createRealIdentityMatrix(Q);
        EssenceMatrix essence = new EssenceMatrix(essenceData);
        essence.setRowMetaData(0, new RowMetaData(5,1));
        essence.setRowMetaData(1, new RowMetaData(5,1));
        essence.setRowMetaData(2, new RowMetaData(5,1));
        essence.setRowMetaData(3, new RowMetaData(5,1));
        if (hasRandom)
        {
            // add a random column with specified mean/variance
            ColumnMetaData randColMD = new ColumnMetaData();
            randColMD.setPredictorType(PredictorType.RANDOM);
            randColMD.setMean(MEAN);
            randColMD.setVariance(VARIANCE);
            essence.setColumnMetaData(3, randColMD);
            
            // build sigma G matrix
            double[][] sigmaG = {{VARIANCE,1,1},{1,VARIANCE,1},{1,1,VARIANCE}};
            params.setSigmaGaussianRandom(new Array2DRowRealMatrix(sigmaG));
            
            // build sigma Y matrix
            double rho = 0.4;
            double [][] sigmaY = {{1,rho,rho},{rho,1,rho},{rho,rho,1}};
            params.setSigmaOutcome(new Array2DRowRealMatrix(sigmaY));

            // build sigma YG
            double rhoYG = 0.8;
            double [][] sigmaYG = {{rhoYG,rho,rho},{rho,rhoYG,rho},{rho,rho,rhoYG}};
            params.setSigmaOutcomeGaussianRandom(new Array2DRowRealMatrix(sigmaYG));
        }
        else
        {
            // build sigma matrix
            double rho = 0.4;
            double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}};
            params.setSigmaError(new Array2DRowRealMatrix(sigma));
        }
        params.setDesignEssence(essence);
        // build beta matrix
        double [][] beta = {{1,0,0},{0,0,0},{0,0,0},{0,0,0}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0,0},{0,0},{0,0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));

        // build between subject contrast
        double [][] between = {{1,-1,0,0},{1,0,-1,0},{1,0,0,-1}};
        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));

        // build within subject contrast
        double [][] within = {{1,1},{-1,0},{0,-1}};
        params.setWithinSubjectContrast(new Array2DRowRealMatrix(within));

        return params;     
    }   
}
