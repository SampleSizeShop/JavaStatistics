package edu.cudenver.bios.powersamplesize;

import jsc.distributions.FishersF;
import jsc.distributions.NoncentralFishersF;
import jsc.distributions.Normal;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.glmm.GLMMTest;
import edu.cudenver.bios.powersamplesize.glmm.GLMMTestFactory;
import edu.cudenver.bios.powersamplesize.glmm.NonCentralityDistribution;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;

/**
 * Calculate power for the general linear multivariate model
 * <p>
 * Power formulas based on the following publication:
 * <p>
 * Muller, K. E., LaVange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
 * Power calculations for general linear multivariate models including repeated 
 * measures applications. Journal of the American Statistical Association, 87: 1209-1226.
 * 
 * @author Sarah Kreidler
 *
 */
public class PowerGLMM implements Power
{       
    /**
     * Calculate power for the general linear multivariate model based on
     * the input matrices.
     * 
     * For the mutivariate case, the following matrices must be specified
     * <ul>
     * <li>Design matrix - essence matrix for the model design
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
    public double getCalculatedPower(PowerSampleSizeParameters params)
            throws IllegalArgumentException
    {
        LinearModelPowerSampleSizeParameters powerParams = (LinearModelPowerSampleSizeParameters) params;
        // make sure all of the matrix inputs are appropriate
        validateParameters(powerParams);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(powerParams);
        
        switch (powerParams.getPowerMethod())
        {
        case CONDITIONAL_POWER:
        	return getConditionalPower(powerParams);
        case QUANTILE_POWER:
        	return getQuantilePower(powerParams);
        case UNCONDITIONAL_POWER:
        	return getUnconditionalPower(powerParams);
        default:
        	return getConditionalPower(powerParams);
        }
    }

    /**
     * Simulate power for the general linear multivariate model based on
     * the input matrices.
     * 
     * For the mutivariate case, the following matrices must be specified
     * <ul>
     * <li>Design matrix - essence matrix for the model design
     * <li>Beta matrix - estimated regression coefficients matrix
     * <li>Sigma matrix - estimated errors matrix
     * <li>Theta0 matrix - estimated null hypothesis matrix 
     * <li>Between subject contrast matrix - defines comparisons between subjects
     * <li>Within subject contrast matrix - defines comparisons within subjects
     * (may be left null for univariate special case)
     * </ul>
     * 
     * @param params Container for input matrices
     * @param interations number of simulation samples/iterations
     * @return simulated power value
     */
    @Override
    public double getSimulatedPower(PowerSampleSizeParameters params, int iterations)
            throws IllegalArgumentException
    {
        LinearModelPowerSampleSizeParameters powerParams = (LinearModelPowerSampleSizeParameters) params;
        validateParameters(powerParams); //  check inputs
     
        // precalculate any computationally expensive matrices, 
        // update the parameters as needed - used for random covariates
        initialize(powerParams);
        
        // get total observations, N, and rank of design matrix
        RealMatrix X = powerParams.getDesign();
        int N = X.getRowDimension();
        int rankX = new SingularValueDecompositionImpl(X).getRank();
                
        // create a normal distribution for generating random errors
        Normal normalDist = new Normal();       
        
        // run the simulations
        int rejectionCount = 0;
        for(int i = 0; i < iterations; i++)
        {
            RealMatrix error = 
                simulateError(normalDist, powerParams.getDesign().getRowDimension(),
                        powerParams.getBeta().getColumnDimension(),
                        powerParams.getSigmaError());         

            // calculate simulated Y based on Y = X beta + error
            RealMatrix Ysim = (X.multiply(powerParams.getBeta())).add(error);
            // calculate beta-Hat
            RealMatrix XtX = X.transpose().multiply(X);
            RealMatrix XtXInverse = new LUDecompositionImpl(XtX).getSolver().getInverse();
            RealMatrix betaHat = XtXInverse.multiply(X.transpose()).multiply(Ysim);
            // calculate Y-hat
            RealMatrix YHat = (X.multiply(betaHat));
            // calculate sigma-Hat
            RealMatrix Ydiff = Ysim.subtract(YHat);
            RealMatrix sigmaHat = (Ydiff.transpose().multiply(Ydiff)).scalarMultiply(((double) 1/(double)(N - rankX)));    
            
            // create a copy of the input parameters, resetting the values of betaHat, sigmaHat,
            // based on the simulated errors
            LinearModelPowerSampleSizeParameters simulatedParams = new LinearModelPowerSampleSizeParameters(powerParams);
            simulatedParams.setBeta(betaHat);
            simulatedParams.setSigmaError(sigmaHat);
            
            // calculate the observed F for the simulation
            GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(simulatedParams);
            double fobs = glmmTest.getObservedF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);

            // get the p-value from a central F distribution
            double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
            double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.DATA_ANALYSIS_NULL);
            double pvalue = FishersF.upperTailProb(fobs, ndf, ddf);

            // check if we reject the null hypothesis
            if (pvalue <= params.getAlpha())
                rejectionCount++;
        }

        return ((double) rejectionCount) / ((double) iterations);
    }

    /**
     * Simulate the error matrix in the Y = X * beta + e
     */
    private RealMatrix simulateError(Normal normalDist, int rows, int columns, RealMatrix sigma)
    throws IllegalArgumentException
    {        
        // build a matrix of random values from a standard normal
        // the number of rows = #subjects (rows) in the full design matrix
        // the number of columns = #outcome variables (i.e. columns in beta)
        Array2DRowRealMatrix randomNormals = new Array2DRowRealMatrix(rows, columns);
        for(int rowIndex = 0; rowIndex < rows; rowIndex++)
        {
            for(int columnIndex = 0; columnIndex < columns; columnIndex++)
            {
                randomNormals.setEntry(rowIndex, columnIndex, normalDist.random()); 
            }
        }
        
        // take the square root of the sigma matrix via cholesky decomposition
        try
        {            
            RealMatrix sqrtMatrix = new CholeskyDecompositionImpl(sigma).getL();
            return randomNormals.multiply(sqrtMatrix); 
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
    }
    
    /**
     * 
     */
    private void validateParameters(LinearModelPowerSampleSizeParameters params)
    throws IllegalArgumentException
    {         
        // convenience variables
        RealMatrix beta = params.getBeta();
        RealMatrix theta0 = params.getTheta();
        RealMatrix X = params.getDesign();
        EssenceMatrix XEssence = params.getDesignEssence();
        int numRandom = (XEssence != null ? XEssence.getRandomPredictorCount() : 0);
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        
        // validate alpha level, 0 < alpha < 1
        if (params.getAlpha() <= 0 || params.getAlpha() >= 1)
            throw new IllegalArgumentException("Alpha must be between 0 and 1.  Invalid value: " + params.getAlpha());
        // make sure sample size > 1
        if (params.getSampleSize() <= 1)
            throw new IllegalArgumentException("Sample size must be greater than 1.  Invalid value: " + params.getSampleSize());
        // only allow at most 1 random predictor
        // TODO: handle multiple random predictors
        if (numRandom > 1)
            throw new IllegalArgumentException("Two many random predictors - at most 1 is allowed");

        // make sure all required matrices have been specified
        // note, we don't check U (within subject contrast), since it may be null in univariate cases
        if (beta == null) 
            throw new IllegalArgumentException("No beta (regression coefficients) matrix specified");
        if (X == null)
            throw new IllegalArgumentException("No design matrix specified");
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
        
        // different variance/covariance matrices are specified depending on the presence
        // of random covariate
        if (numRandom == 0)
        {
            if (sigmaE == null)
                throw new IllegalArgumentException("No sigma (error) matrix specified");
            if (!sigmaE.isSquare())
                throw new IllegalArgumentException("Sigma error matrix must be square");
            if (U.getRowDimension() != sigmaE.getRowDimension())
                throw new IllegalArgumentException("Within subject contrast does not conform with sigma matrix");
        }
        else if (numRandom == 1)
        {
            // make sure the test statistic is either HLT or UNIREP if there is a random
            // covariate (results not published for Wilk's Lambda or Pillai-Bartlett 
            if (params.getTestStatistic() != TestStatistic.HOTELLING_LAWLEY_TRACE &&
                    params.getTestStatistic() != TestStatistic.UNIREP)
                throw new IllegalArgumentException("With a random covariate, only Hotelling-Lawley and Unirep test statistics are supported");
            
            if (sigmaG == null)
                throw new IllegalArgumentException("No variance/covariance matrix specified for gaussian predictors");
            if (sigmaY == null)
                throw new IllegalArgumentException("No variance/covariance matrix specified for response variables");
            if (sigmaYG == null)
                throw new IllegalArgumentException("No outcome / gaussian predictor covariance matrix specified");
            
            // check conformance
            if (U.getRowDimension() != sigmaY.getRowDimension())
                throw new IllegalArgumentException("Within subject contrast does not conform with sigma matrix");
            if (sigmaG.getRowDimension() != sigmaYG.getColumnDimension())
                throw new IllegalArgumentException("Outcome / Gaussian predictor covariance matrix does not conform with variance matrix for the gaussian predictor");
            if (!sigmaY.isSquare())
                throw new IllegalArgumentException("Variance/covariance matrix for response variables must be square");
            if (!sigmaG.isSquare())
                throw new IllegalArgumentException("Variance/covariance matrix for gaussian predictors must be square");
        }
        
        // check dimensionality 
        if (C.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Between subject contrast does not conform with beta matrix");
        if (beta.getColumnDimension() != U.getRowDimension())
            throw new IllegalArgumentException("Within subject contrast does not conform with beta matrix");
        if (X.getColumnDimension() != beta.getRowDimension())
            throw new IllegalArgumentException("Design matrix does not conform with beta matrix");
        if (C.getRowDimension() > C.getColumnDimension())
            throw new IllegalArgumentException("Number of rows in between subject contrast must be less than or equal to the number of columns");
        if (U.getColumnDimension() > U.getRowDimension())
            throw new IllegalArgumentException("Number of columns in within subject contrast must be less than or equal to the number of rows");
        if (theta0.getRowDimension() != C.getRowDimension())
            throw new IllegalArgumentException("Number of rows in theta null must equal number of rows in between subject contrast");

        // check rank of the design matrix
        int rankX = new SingularValueDecompositionImpl(X).getRank();
        if (rankX != X.getColumnDimension())
            throw new IllegalArgumentException("Design matrix is not full rank");

        // make sure design matrix is symmetric and positive definite
        // TODO: how to check this?
    }
    
    private void initialize(LinearModelPowerSampleSizeParameters params)
    {        
        // update the sigma error if we have a baseline covariate
        EssenceMatrix XEssence = params.getDesignEssence();
        int numRandom = (XEssence != null ? XEssence.getRandomPredictorCount() : 0);
        if (numRandom == 1)
        {
            RealMatrix sigmaG = params.getSigmaGaussianRandom();
            RealMatrix sigmaY = params.getSigmaOutcome();
            RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
            
            // set the sigma error matrix to [sigmaY - sigmaYG * sigmaG-1 * sigmaGY] 
            RealMatrix sigmaGY = sigmaYG.transpose();
            RealMatrix sigmaGInverse = new LUDecompositionImpl(sigmaG).getSolver().getInverse();
            params.setSigmaError(sigmaY.subtract(sigmaYG.multiply(sigmaGInverse.multiply(sigmaGY))));
            
            // calculate the betaG matrix and fill in the placeholder row for the random predictor
            RealMatrix beta = params.getBeta();
            // first, find the random predictor column index
            // TODO: maybe a more convenient function on EssenceMatrix class?
            int randomCol = -1;
            for (int col = 0; col < XEssence.getColumnDimension(); col++)
            {
                ColumnMetaData colMD = XEssence.getColumnMetaData(col);
                if (colMD.getPredictorType() == PredictorType.RANDOM)
                {
                    randomCol = col;
                    break;
                }
            }
            RealMatrix betaG = sigmaGInverse.multiply(sigmaGY);
            for (int col = 0; col < betaG.getColumnDimension(); col++)
            {
                beta.setEntry(randomCol, col, betaG.getEntry(0, col));
            }
        }
    }
    
    
    private double getConditionalPower(LinearModelPowerSampleSizeParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic 
        // under the null hypothesis
        double nonCentralityParam = glmmTest.getNonCentrality(GLMMTest.DistributionType.POWER_ALTERNATIVE);

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }
    
    private double getUnconditionalPower(LinearModelPowerSampleSizeParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getAlpha());
        
        // determine unconditional power by integrating over all possible values 
        // of the non-centrality parameter
        // TODO: write me :-)
        
        return 0;
        
    }
    
    private double getQuantilePower(LinearModelPowerSampleSizeParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        NonCentralityDistribution nonCentralityDist = 
            new NonCentralityDistribution(params, false);
        double nonCentralityParam = nonCentralityDist.inverseCDF(params.getQuantile());
        
        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        double altDdf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }

}
