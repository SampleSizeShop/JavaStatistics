package edu.cudenver.bios.powersamplesize;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;

/**
 * Calculate sample size for the general linear multivariate model
 * 
 * @author Sarah Kreidler
 *
 */
public class SampleSizeGLMM implements SampleSize
{
    private static final int STARTING_SAMPLE_SIZE = 1000;
    private static final int MIN_SAMPLE_SIZE =  2; // need df > 0
    
    // function to be used with apache's built-in bisection solver
    private class SampleSizeFunction implements UnivariateRealFunction
    {
        LinearModelPowerSampleSizeParameters params;
        PowerGLMM powerGLMM;
        
        public SampleSizeFunction(LinearModelPowerSampleSizeParameters params)
        {
            this.params = params;
            this.powerGLMM = new PowerGLMM();
        }
        
        public double value(double n)
        {
            try
            {
                LinearModelPowerSampleSizeParameters tmpParams = 
                    new LinearModelPowerSampleSizeParameters(params);
                RealMatrix design = tmpParams.getDesignEssence().getFullDesignMatrix((int) n);
                tmpParams.setDesign(design);
                tmpParams.setSampleSize((int) n);
                double calculatedPower = powerGLMM.getCalculatedPower(tmpParams);
                return tmpParams.getPower() - calculatedPower;
            }
            catch (Exception e)
            {   
                // we can't throw an exception here since the UnivariateRealFunction interface does
                // not allow it.  So we return a large number to prevent BisectionSolver from using
                // the n which caused to exception as a root
                return STARTING_SAMPLE_SIZE;  
            }
        }
    }
    
    /**
     * Calculate power for the general linear multivariate model based on
     * the input matrices.
     * 
     * For the mutivariate case, the following matrices must be specified
     * <ul>
     * <li>Essence design matrix - essence matrix for the model design
     * <li>Beta matrix - estimated regression coefficients matrix
     * <li>Sigma matrix - estimated errors matrix
     * <li>Theta0 matrix - estimated null hypothesis matrix 
     * <li>Between subject contrast matrix - defines comparisons between subjects
     * <li>Within subject contrast matrix - defines comparisons within subjects
     * (may be left null for univariate special case)
     * </ul>
     * 
     * @param params Container for input matrices
     * @return power
     */
    @Override
    public int getSampleSize(PowerSampleSizeParameters params)
            throws IllegalArgumentException
    {
        LinearModelPowerSampleSizeParameters ssParams =
            (LinearModelPowerSampleSizeParameters) params;
        validateParameters(ssParams);

        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        SampleSizeFunction sampleSizeFunc = new SampleSizeFunction(ssParams);
        
        try
        {
            int upperBound = getUpperBound(ssParams);
            int minSampleSize = ssParams.getDesignEssence().getMinimumSampleSize();
            int adjSize = (int) Math.ceil(solver.solve(sampleSizeFunc, minSampleSize, upperBound));
            // increment the calculated sample size to get a value
            // that is a multiple of the minimum sample size
            int remainder = adjSize % minSampleSize;
            if (remainder != 0) adjSize += minSampleSize - remainder;
            return adjSize;
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to calculate sample size: " + e.getMessage());
        }
    }

    private int getUpperBound(LinearModelPowerSampleSizeParameters params)
    {
        PowerGLMM calc = new PowerGLMM();
        LinearModelPowerSampleSizeParameters tmpParams = 
            new LinearModelPowerSampleSizeParameters(params);
        double desiredPower = params.getPower();
        int upperBound = STARTING_SAMPLE_SIZE;
        
        for(double currentPower = 0.0; currentPower < desiredPower; upperBound *= 2)
        {
            RealMatrix design = params.getDesignEssence().getFullDesignMatrix(upperBound);
            tmpParams.setDesign(design);
            tmpParams.setSampleSize(upperBound);
            currentPower = calc.getCalculatedPower(tmpParams);
        }
        return upperBound;
    }
    
    /**
     * Validate the incoming linear model parameters
     * 
     * @param params
     * @throws IllegalArgumentException
     */
    private void validateParameters(LinearModelPowerSampleSizeParameters params)
    throws IllegalArgumentException
    {
        // convenience variables
        RealMatrix beta = params.getBeta();
        RealMatrix sigma = params.getSigmaError();
        RealMatrix theta0 = params.getTheta();
        EssenceMatrix essenceX = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        
        // validate alpha level, 0 < alpha < 1
        if (params.getAlpha() <= 0 || params.getAlpha() >= 1)
            throw new IllegalArgumentException("Alpha must be between 0 and 1.  Invalid value: " + params.getAlpha());
        // make sure power > alpha and < 1
        if (params.getPower() <= params.getAlpha() || params.getPower() >= 1)
            throw new IllegalArgumentException("Power must be between alpha and 1.  Invalid value: " + params.getPower());
        
        // make sure all required matrices have been specified
        // note, we don't check U (within subject contrast), since it may be null in univariate cases
        if (beta == null) 
            throw new IllegalArgumentException("No beta (regression coefficients) matrix specified");
        if (essenceX == null)
            throw new IllegalArgumentException("No design essence matrix specified");
        if (sigma== null)
            throw new IllegalArgumentException("No sigma (error) matrix specified");
        if (C == null)
            throw new IllegalArgumentException("No between subject contrast (C) matrix specified");
        if (theta0 == null)
            throw new IllegalArgumentException("No theta_null (null hypothesis) matrix specified");
        // create a default U if not specified
        if (U == null)
        {
            U = MatrixUtils.createRealIdentityMatrix(beta.getColumnDimension());
            params.setWithinSubjectContrast(U);
        }
        
        // check dimensionality 
        if (U.getRowDimension() != sigma.getRowDimension())
            throw new IllegalArgumentException("Within subject contrast does not conform with sigma matrix");
        if (C.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Between subject contrast does not conform with beta matrix");
        if (beta.getColumnDimension() != U.getRowDimension())
            throw new IllegalArgumentException("Within subject contrast does not conform with beta matrix");
        if (essenceX.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Design matrix does not conform with beta matrix");
        if (C.getRowDimension() > C.getColumnDimension())
            throw new IllegalArgumentException("Number of rows in between subject contrast must be less than or equal to the number of columns");
        if (U.getColumnDimension() > U.getRowDimension())
            throw new IllegalArgumentException("Number of columns in within subject contrast must be less than or equal to the number of rows");
        if (theta0.getRowDimension() != C.getRowDimension())
            throw new IllegalArgumentException("Number of rows in theta null must equal number of rows in between subject contrast");

        // make sure at least one of the entries in the beta matrix is different than the
        // others
        if (beta.getColumnDimension() <= 1 && beta.getRowDimension() <= 1)
            throw new IllegalArgumentException("Must have at least two groups for sample size calculation");
        double betaValue = beta.getEntry(0, 0);
        boolean foundDiff = false;
        for(int row = 0; row < beta.getRowDimension() && !foundDiff; row++)
        {
            for(int col = 0; col < beta.getColumnDimension() & !foundDiff; col++)
            {
                if (betaValue != beta.getEntry(row, col)) foundDiff = true;
            }
        }
        if (foundDiff == false) 
            throw new IllegalArgumentException("Must specify at least one different mean for sample size calculation");
    }

}
