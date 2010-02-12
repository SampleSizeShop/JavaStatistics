package edu.cudenver.bios.powersamplesize.test;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.PowerGLMM;
import edu.cudenver.bios.powersamplesize.PowerOneSampleStudentsT;
import edu.cudenver.bios.powersamplesize.SampleSizeGLMM;
import edu.cudenver.bios.powersamplesize.SampleSizeOneSampleStudentsT;
import edu.cudenver.bios.powersamplesize.graphics.PowerCurveBuilder;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.SimplePowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;
import junit.framework.TestCase;

public class TestPowerCurveBuilder extends TestCase
{
    private static final double ALPHA = 0.05;
    private static final int SAMPLE_SIZE = 10;
    private static final String OUTPUT_DIR = "C:\\Documents and Settings\\kreidles\\Desktop\\";
    
    private void testPowerCurveOneSampleStudentsTbyN()
    {
        SimplePowerSampleSizeParameters params = 
            new SimplePowerSampleSizeParameters();
        params.setMu0(0);
        params.setMuA(2);
        params.setSigma(1.5);
        params.setAlpha(0.05);
        params.setSampleSize(10);
        PowerCurveBuilder curveBuilder = 
            new PowerCurveBuilder(new PowerOneSampleStudentsT(), 
                    new SampleSizeOneSampleStudentsT());
        curveBuilder.setLegend(true);
        curveBuilder.setTitle("PowerCurve Unit Test");
        curveBuilder.setXaxisLabel("Total Sample Size");
        curveBuilder.setYaxisLabel("Power");
        JFreeChart chart = curveBuilder.getPowerCurve(params);
        
        // write the chart as a jpeg image stream 
        try
        {
            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "studentTbyN.jpeg"), chart, 500, 300);
        }
        catch (IOException e)
        {
            System.err.println("write failed: " + e.getMessage());
        }
    }
    
    private void testPowerCurveOneSampleStudentsTbyDelta()
    {
        SimplePowerSampleSizeParameters params = 
            new SimplePowerSampleSizeParameters();
        params.setMu0(0);
        params.setMuA(2);
        params.setSigma(1.5);
        params.setAlpha(0.05);
        params.setSampleSize(10);
        PowerCurveBuilder curveBuilder = 
            new PowerCurveBuilder(new PowerOneSampleStudentsT(), 
                    new SampleSizeOneSampleStudentsT());
        curveBuilder.setBySampleSize(false);
        curveBuilder.setLegend(true);
        curveBuilder.setMinimumMeanDifference(0);
        curveBuilder.setMeanDifferenceIncrement(0.05);
        curveBuilder.setTitle("PowerCurve Unit Test");
        curveBuilder.setXaxisLabel("Mean Difference");
        curveBuilder.setYaxisLabel("Power");
        JFreeChart chart = curveBuilder.getPowerCurve(params);
        
        // write the chart as a jpeg image stream 
        try
        {
            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "studentTbydelta.jpeg"), chart, 500, 300);
        }
        catch (IOException e)
        {
            System.err.println("write failed: " + e.getMessage());
        }
    }
    
    private void testPowerCurveGLMMUnivariateByN()
    {
        LinearModelPowerSampleSizeParameters params = buildValidUnivariateInputs();
        params.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        PowerCurveBuilder curveBuilder = 
            new PowerCurveBuilder(new PowerGLMM(), new SampleSizeGLMM());
        curveBuilder.setMinimumSampleSize(2*params.getDesignEssence().getMinimumSampleSize());
        curveBuilder.setSampleSizeIncrement(params.getDesignEssence().getMinimumSampleSize());
        curveBuilder.setLegend(true);
        curveBuilder.setTitle("PowerCurve Unit Test: GLMM Univariate Case");
        curveBuilder.setXaxisLabel("Total Sample Size");
        curveBuilder.setYaxisLabel("Power");
        JFreeChart chart = curveBuilder.getPowerCurve(params);
        
        // write the chart as a jpeg image stream 
        try
        {
            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "univariateN.jpeg"), chart, 500, 300);
        }
        catch (IOException e)
        {
            System.err.println("write failed: " + e.getMessage());
        }
    }
    
    private void testPowerCurveGLMMUnivariateByDelta()
    {
        LinearModelPowerSampleSizeParameters params = buildValidUnivariateInputs();
        params.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        PowerCurveBuilder curveBuilder = 
            new PowerCurveBuilder(new PowerGLMM(), new SampleSizeGLMM());
        curveBuilder.setBySampleSize(false);
        curveBuilder.setLegend(true);
        curveBuilder.setTitle("PowerCurve Unit Test: GLMM Univariate Case");
        curveBuilder.setXaxisLabel("Mean Difference");
        curveBuilder.setYaxisLabel("Power");
        JFreeChart chart = curveBuilder.getPowerCurve(params);
        
        // write the chart as a jpeg image stream 
        try
        {
            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "univariateDelta.jpeg"), chart, 500, 300);
        }
        catch (IOException e)
        {
            System.err.println("write failed: " + e.getMessage());
        }
    }
    
    public void testPowerCurveGLMMMultivariateByN()
    {
        LinearModelPowerSampleSizeParameters params = buildValidMultivariateFixedInputs();
        params.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
        PowerCurveBuilder curveBuilder = 
            new PowerCurveBuilder(new PowerGLMM(), new SampleSizeGLMM());
        curveBuilder.setBySampleSize(true);
        curveBuilder.setLegend(true);
        curveBuilder.setTitle("PowerCurve Unit Test: GLMM Multivariate Case");
        curveBuilder.setXaxisLabel("N");
        curveBuilder.setYaxisLabel("Power");
        JFreeChart chart = curveBuilder.getPowerCurve(params);
        
        // write the chart as a jpeg image stream 
        try
        {
            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "multivariateN.jpeg"), chart, 500, 300);
        }
        catch (IOException e)
        {
            System.err.println("write failed: " + e.getMessage());
        }
    }
    
//    public void testPowerCurveGLMMMultivariateByDelta()
//    {
//        LinearModelPowerSampleSizeParameters params = buildValidUnivariateInputs();
//        params.setTestStatistic(TestStatistic.HOTELLING_LAWLEY_TRACE);
//        PowerCurveBuilder curveBuilder = 
//            new PowerCurveBuilder(new PowerGLMM(), new SampleSizeGLMM());
//        curveBuilder.setBySampleSize(false);
//        curveBuilder.setLegend(true);
//        curveBuilder.setTitle("PowerCurve Unit Test: GLMM Univariate Case");
//        curveBuilder.setXaxisLabel("Mean Difference");
//        curveBuilder.setYaxisLabel("Power");
//        JFreeChart chart = curveBuilder.getPowerCurve(params);
//        
//        // write the chart as a jpeg image stream 
//        try
//        {
//            ChartUtilities.saveChartAsJPEG(new File(OUTPUT_DIR + "univariateDelta.jpeg"), chart, 500, 300);
//        }
//        catch (IOException e)
//        {
//            System.err.println("write failed: " + e.getMessage());
//        }
//    }
    
    private LinearModelPowerSampleSizeParameters buildValidUnivariateInputs()
    {
        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
        params.setAlpha(ALPHA);
        
        // build beta matrix
        double [][] beta = {{0},{1.05}};
        params.setBeta(new Array2DRowRealMatrix(beta));
        // build theta null matrix
        double [][] theta0 = {{0}};
        params.setTheta(new Array2DRowRealMatrix(theta0));
        // build sigma matrix
        double [][] sigma = {{2.05}};
        params.setSigmaError(new Array2DRowRealMatrix(sigma));
        // build design matrix

        double[][] essenceData = {{1,0},{0,1}};
        EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
        essenceMatrix.setRowMetaData(0, new RowMetaData(10,1));
        essenceMatrix.setRowMetaData(1, new RowMetaData(10,1));
//        ColumnMetaData cmd = new ColumnMetaData();
//        cmd.setMean(5);
//        cmd.setVariance(3);
//        cmd.setPredictorType(PredictorType.RANDOM);
//        essenceMatrix.setColumnMetaData(1, cmd);
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
}
