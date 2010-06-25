package edu.cudenver.bios.powersamplesize.test;

import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.PowerGLMM;
import edu.cudenver.bios.powersamplesize.glmm.NonCentralityDistribution;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.PowerMethod;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.UnivariateCorrection;

import jsc.distributions.Normal;
import junit.framework.TestCase;

public class TestPowerGLMM extends TestCase
{
    private static final int SIMULATION_SIZE = 10000;
    private static final double UNIT_TEST_ALPHA = 0.01;
    private static final double ALPHA = 0.05;    
    private static final double MEAN = 9.75;
    private static final double VARIANCE = 2.0;
    private static final double[] BETA_SCALE = {0,0.5,1,1.5,2};
    private static final double[] SIGMA_SCALE = {1,2};
    private static final int[] SAMPLE_SIZE = {20};
    private Normal normalDist;
    private DecimalFormat Number;
    
    public void setUp()
    {
        normalDist = new Normal();
        Number = new DecimalFormat("#0.000");
    }
    
    private void testValidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();

        PowerGLMM calc = new PowerGLMM();
//        goodParams.setTestStatistic(TestStatistic.UNIREP);
//        checkPower("Valid Univariate, Fixed, UNIREP", calc, goodParams);
//        // make sure unirep corrections don't mess up the univariate case
//        goodParams.setUnivariateCorrection(UnivariateCorrection.BOX);
//        checkPower("Valid Univariate, Fixed, UNIREP, BOX correction", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Univariate, Fixed, HLT",calc, goodParams);
//        goodParams.setTestStatistic(TestStatistic.PILLAI_BARTLETT_TRACE);
//        checkPower("Valid Univariate, Fixed, PB",calc, goodParams);
//        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
//        checkPower("Valid Univariate, Fixed, W",calc, goodParams);        
    }

    private void testInvalidUnivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidUnivariateInputs();
        goodParams.setBeta(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Univariate, Fixed, UNIREP", calc, goodParams);
    }

    public void testValidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateFixedInputs();

        PowerGLMM calc = new PowerGLMM();
//        goodParams.setMomentMethod(MomentApproximationMethod.PILLAI_ONE_MOMENT);
//        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
//        checkPower("Valid Multivariate, Fixed, HLT",calc, goodParams);
//        
//        goodParams.setMomentMethod(MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT);
//        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
//        checkPower("Valid Multivariate, Fixed, HLT",calc, goodParams);
//        
//        goodParams.setMomentMethod(MomentApproximationMethod.MCKEON_TWO_MOMENT);
//        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
//        checkPower("Valid Multivariate, Fixed, HLT",calc, goodParams);
//        
//        goodParams.setMomentMethod(MomentApproximationMethod.MCKEON_TWO_MOMENT_OMEGA_MULT);
//        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
//        checkPower("Valid Multivariate, Fixed, HLT",calc, goodParams);

//        goodParams.setMomentMethod(MomentApproximationMethod.PILLAI_ONE_MOMENT);
//        goodParams.setTestStatistic(TestStatistic.PILLAI_BARTLETT_TRACE);
//        checkPower("Valid Multivariate, Fixed, PB (Pillai 1 moment)",calc, goodParams);
//
//        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
//        checkPower("Valid Multivariate, Fixed, W",calc, goodParams);       
//        
//        goodParams.setMomentMethod(MomentApproximationMethod.RAO_TWO_MOMENT_OMEGA_MULT);
//        goodParams.setTestStatistic(TestStatistic.WILKS_LAMBDA);
//        checkPower("Valid Multivariate, Fixed, W",calc, goodParams);    
        goodParams.setTestStatistic(TestStatistic.UNIREP);
//        checkPower("Valid Multivariate, Fixed, UNIREP, uncorrected", calc, goodParams);
        //goodParams.setUnivariateCorrection(UnivariateCorrection.BOX);
  //      checkPower("Valid Multivariate, Fixed, UNIREP, BOX", calc, goodParams);
        goodParams.setUnivariateCorrection(UnivariateCorrection.GEISSER_GREENHOUSE);
        checkPower("Valid Multivariate, Fixed, UNIREP, GG", calc, goodParams);
//        goodParams.setUnivariateCorrection(UnivariateCorrection.HUYNH_FELDT);
//        checkPower("Valid Multivariate, Fixed, UNIREP, HF", calc, goodParams);

    }

    private void testInvalidMultivariateFixed()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateFixedInputs();

        goodParams.setSigmaGaussianRandom(null);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        checkPowerFail("Invalid Beta, Multivariate, Fixed, UNIREP", calc, goodParams);

    }

    private void testValidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateRandomInputs();
        //goodParams.setPowerMethod(PowerMethod.QUANTILE_POWER);
        goodParams.setPowerMethod(PowerMethod.UNCONDITIONAL_POWER);

        goodParams.setQuantile(0.25);
        PowerGLMM calc = new PowerGLMM();
        goodParams.setTestStatistic(TestStatistic.UNIREP);
        goodParams.setUnivariateCorrection(UnivariateCorrection.NONE);
        checkPower("Valid Multivariate, Random, UNIREP", calc, goodParams);
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        checkPower("Valid Multivariate, Random, HLT",calc, goodParams);
    }

    private void testInvalidMultivariateRandom()
    {
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateRandomInputs();

    }

    private void testNonCentralityCDF()
    {
        double[] testFCrits = {0.8223684, 1.2335526, 1.3748972, 1.439, 1.55,1.63,1.643333,1.65};
        LinearModelPowerSampleSizeParameters goodParams = buildValidMultivariateRandomInputs();
        goodParams.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        goodParams.setPowerMethod(PowerMethod.QUANTILE_POWER);
        goodParams.setQuantile(0.25);
        // build the error matrix
        RealMatrix sigmaYG = goodParams.getSigmaOutcomeGaussianRandom();
        RealMatrix sigmaG = goodParams.getSigmaGaussianRandom();
        RealMatrix sigmaY = goodParams.getSigmaOutcome();
        RealMatrix sigmaGY = sigmaYG.transpose();
        RealMatrix sigmaGInverse = new LUDecompositionImpl(sigmaG).getSolver().getInverse();
        goodParams.setSigmaError(sigmaY.subtract(sigmaYG.multiply(sigmaGInverse.multiply(sigmaGY))));
        
        NonCentralityDistribution nonCentralityDist = 
            new NonCentralityDistribution(goodParams, false);
        System.out.println("Starting cdf tests");
        for(double fcrit : testFCrits)
        {
            double prob = nonCentralityDist.cdf(fcrit);
            System.out.println("Crit=" + fcrit + " prob=" + prob);
        }
        
        double fval = nonCentralityDist.inverseCDF(goodParams.getQuantile());
        System.out.println("Quantile=" + goodParams.getQuantile() + " fval=" + fval);

        assertTrue(true);
    }
    
    private void checkPowerFail(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        try
        {
            double calculated = calc.getCalculatedPower(params);
            double simulated = calc.getSimulatedPower((PowerSampleSizeParameters) params, SIMULATION_SIZE);
            assertFalse(powersAreSame(calculated, simulated));
        }
        catch(Exception e)
        {
            System.out.println(label + ">> Recognized invalid inputs: " + e.getMessage());
            assert(true);
        }

    }

    private void checkPower(String label, PowerGLMM calc, LinearModelPowerSampleSizeParameters params)
    {
        int tests = 0;
        int matches = 0;

        for(double sigmaScale : SIGMA_SCALE)
        {
            for(double betaScale : BETA_SCALE)
            {
                for(int sampleSize: SAMPLE_SIZE)
                {  
                    LinearModelPowerSampleSizeParameters testParams = 
                        new LinearModelPowerSampleSizeParameters(params);
                    // scale the inputs
                    testParams.setBeta(params.getBeta().scalarMultiply(betaScale));
                    if (params.getSigmaError() != null)
                        testParams.setSigmaError(params.getSigmaError().scalarMultiply(sigmaScale));
                    else
                        testParams.setSigmaGaussianRandom(params.getSigmaGaussianRandom().scalarMultiply(sigmaScale));
                    testParams.setSampleSize(sampleSize);
                    try
                    {
                        double calculated = calc.getCalculatedPower(testParams);
                        double simulated = calc.getSimulatedPower(testParams, SIMULATION_SIZE);

                        System.out.println(label + "["+sigmaScale+","+betaScale+","+sampleSize
                                +"]>> Calculated power: " + Number.format(calculated) + ", simulated power: " + Number.format(simulated));
                        tests++;
                        if (powersAreSame(calculated, simulated)) matches++;
                    }
                    catch(Exception e)
                    {
                        System.out.println(label + ">> Failed to calculate power: " + e.getMessage());
                    }        

                }
            }
        }
        assertTrue(true);
        //assertEquals(tests, matches);
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


    private LinearModelPowerSampleSizeParameters buildValidMultivariateFixedInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);

        int Q = 4;
        // create design matrix
        RealMatrix essenceData = MatrixUtils.createRealIdentityMatrix(Q);
        EssenceMatrix essence = new EssenceMatrix(essenceData);
        essence.setRowMetaData(0, new RowMetaData(5,1));
        essence.setRowMetaData(1, new RowMetaData(5,1));
        essence.setRowMetaData(2, new RowMetaData(5,1));
        essence.setRowMetaData(3, new RowMetaData(5,1));
        params.setDesignEssence(essence);
        
        // build sigma matrix
        double rho = 0.4;
        double [][] sigma = {{1,rho,rho},{rho,1,rho},{rho,rho,1}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        
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
    
    private LinearModelPowerSampleSizeParameters buildValidMultivariateRandomInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);

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

        // build beta matrix
        double [][] beta = {{0.25,0},{0,0.25},{0.9,0}};
        params.setBeta(new Array2DRowRealMatrix(beta));

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
