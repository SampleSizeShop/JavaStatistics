package edu.cudenver.bios.powersamplesize.test;

import org.apache.commons.math.linear.Array2DRowRealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.RowMetaData;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;

public class ParameterBuilder
{
	   public static LinearModelPowerSampleSizeParameters buildValidUnivariateInputs(double alpha, int nPerGroup, boolean essenceOnly)
	    {
	        LinearModelPowerSampleSizeParameters params = new LinearModelPowerSampleSizeParameters();
	        params.setAlpha(alpha);
	        params.setSampleSize(10);
	        
	        // build beta matrix
	        double [][] beta = {{0},{1.5}};
	        params.setBeta(new Array2DRowRealMatrix(beta));
	        // build theta null matrix
	        double [][] theta0 = {{0}};
	        params.setTheta(new Array2DRowRealMatrix(theta0));
	        // build sigma matrix
	        double [][] sigma = {{2.05}};
	        params.setSigma(new Array2DRowRealMatrix(sigma));
	        // build design matrix
	        if (!essenceOnly)
	        {
	            double [][] design = {
	                    {1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},{1,0},
	                    {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}};
	            params.setDesign(new Array2DRowRealMatrix(design));
	        }
	        else
	        {
	            double[][] essenceData = {{1,0},{0,1}};
	            EssenceMatrix essenceMatrix = new EssenceMatrix(essenceData);
	            essenceMatrix.setRowMetaData(0, new RowMetaData(5,1));
	            essenceMatrix.setRowMetaData(1, new RowMetaData(5,1));
	            ColumnMetaData cmd = new ColumnMetaData();
	            cmd.setMean(5);
	            cmd.setVariance(3);
	            cmd.setPredictorType(PredictorType.RANDOM);
	            essenceMatrix.setColumnMetaData(1, cmd);
	            params.setDesignEssence(essenceMatrix);
	        }
	        // build between subject contrast
	        double [][] between = {{1,-1}};
	        params.setBetweenSubjectContrast(new Array2DRowRealMatrix(between));
	        
	        return params;     
	        
	    }
}
