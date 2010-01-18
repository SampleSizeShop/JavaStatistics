package edu.cudenver.bios.powersamplesize;

import jsc.distributions.FishersF;
import jsc.distributions.NoncentralFishersF;
import jsc.distributions.Normal;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.PowerMethod;

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
        
        // update the parameters as needed - used for random covariates
        updateParameters(powerParams);
        
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
        
        // update parameters as needed - used for random predictors to set sigmaE, betaG
        updateParameters(powerParams);
        
        // calculate numerator degrees of freedom - these don't change with simulation
        RealMatrix C = powerParams.getBetweenSubjectContrast();
        RealMatrix U = powerParams.getWithinSubjectContrast();
        // if multivariate, then numerator df is a*b, for univariate then numerator df just a.
        int ndf = C.getRowDimension() * U.getColumnDimension();
        // denominator df depends on test, but none of the simulated matrices
        int ddf = getDenominatorDF(powerParams);  
        
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
            double fobs = getObservedF(simulatedParams, ndf, ddf);

            // get the p-value from a central F distribution
            FishersF centralFDist = new FishersF(ndf, ddf);
            double pvalue = 1 - centralFDist.cdf(fobs);

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
    
    private void updateParameters(LinearModelPowerSampleSizeParameters params)
    {
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
    
    
    private double getConditionalPower(LinearModelPowerSampleSizeParameters powerParams)
    {
        // calculate degrees of freedom
        RealMatrix C = powerParams.getBetweenSubjectContrast();
        RealMatrix U = powerParams.getWithinSubjectContrast();
        // if multivariate, then numerator df is a*b, for univariate then numerator df just a.
        int ndf = C.getRowDimension() * U.getColumnDimension();
        // denominator df depends on test
        int ddf = getDenominatorDF(powerParams);
        
        // get the approximate critical F value from a central F distribution
        FishersF centralFDist = new FishersF(ndf, ddf);
        double Fcrit = centralFDist.inverseCdf(1 - powerParams.getAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic
        double nonCentralityParam = ndf * getObservedF(powerParams, ndf, ddf);
        
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(ndf, ddf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }
    
    private double getUnconditionalPower(LinearModelPowerSampleSizeParameters params)
    {
    	return 0.0;
    }
    
    private double getQuantilePower(LinearModelPowerSampleSizeParameters params)
    {
    	return 0.0;
    }
    
    /**
     * Calculate the denominator degrees of freedom for the specified test
     * 
     * @param powerParams input matrices
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    private int getDenominatorDF(LinearModelPowerSampleSizeParameters powerParams)
    throws IllegalArgumentException
    {
        RealMatrix X = powerParams.getDesign();
        RealMatrix C = powerParams.getBetweenSubjectContrast();
        RealMatrix U = powerParams.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        int a = C.getRowDimension();
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();
        // N = total number of subjects (rows in design matrix, X)
        int N = X.getRowDimension();
        // r = rank of design matrix, X
        int r = new SingularValueDecompositionImpl(X).getRank();
        // minimum of a and b dimensions
        int s = (a < b) ? a : b;  
        
        int df = -1;
        switch (powerParams.getTestStatistic())
        {
        case WILKS_LAMBDA:
        {
            if (a*a*b*b <= 4)
            {
                df = N - r - b + 1;
            }
            else
            {
                double gDenominator = (a*a + b*b - 5);
                if (gDenominator == 0)
                    throw new IllegalArgumentException("Within and between subject contrasts yielded divide by zero: row of C=" + a + ", cols of U=" + b);
                double g = Math.sqrt((a*a*b*b - 4) / gDenominator);
                df = (int) Math.round(g*((N - r) - (b - a +1)/2) - (a*b - 2)/2);
            }
            break;
        }
        case PILLAU_BARTLETT_TRACE:
        {
            df = s * ((N - r) - b + s);
            break;
        }
        case UNIREP:
            df = b*(N - r);
            break;
        case HOTELLING_LAWLEY_TRACE:
        {
            df = s * ((N - r) - b -1) + 2;
            break;
        }
        default:
            throw new IllegalArgumentException("Unknown test statistic");
        }
        if (df <= 0)
            throw new IllegalArgumentException("Non-postive denominator degrees of freedom: " + df);
        
        return df;
    }
    
    /**
     * Compute the observed F value under the alternative hypothesis
     * 
     * @param params matrix inputs
     * @param denominator degrees of freedom
     * @return non-centrality parameter for the specific test statistic
     */
    private double getObservedF(LinearModelPowerSampleSizeParameters powerParams, int ndf, int ddf)
    throws IllegalArgumentException
    {
        
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(powerParams);
        RealMatrix errorSumOfSquares = getErrorSumOfSquares(powerParams);
        
        RealMatrix C = powerParams.getBetweenSubjectContrast();
        RealMatrix U = powerParams.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        int a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        int b = U.getColumnDimension();
       // minimum of a and b dimensions
        int s = (a < b) ? a : b;  
        
        double association = 0.0;
        
        switch (powerParams.getTestStatistic())
        {
        case WILKS_LAMBDA:
            double W = getWilksLambda(hypothesisSumOfSquares, errorSumOfSquares);
            if (a*a*b*b <= 4) 
            {
                association = 1 - W;
            }
            else
            {
                double g = Math.sqrt((a*a*b*b - 4) / (a*a + b*b - 5));
                association = 1 - Math.pow(W, 1/g);
            }
            break;
        case PILLAU_BARTLETT_TRACE:
            double PB = getPillauBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
            association = PB / s;
            break;
        case UNIREP:
            double REP = getUnirep(hypothesisSumOfSquares, errorSumOfSquares);
            association = REP / (1 + REP);
            break;
        case HOTELLING_LAWLEY_TRACE:
            double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
            association = (HLT/s) / (1 + (HLT/s));
            break;
        default:
            throw new IllegalArgumentException("Unknown statistic, cannot compute non-centrality parameter.");
        }
        
        double fobs = ((association) / ndf)/ ((1 - association) / ddf);
        return fobs;
    }
    
    /**
     * Compute a Wilks Lamba statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getWilksLambda(RealMatrix H, RealMatrix E)
    throws InvalidMatrixException
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Wilks Lambda: hypothesis and error matrices must be square and same dimensions");
        
        RealMatrix T = H.add(E);
        RealMatrix inverseT = new LUDecompositionImpl(T).getSolver().getInverse();

        RealMatrix EinverseT = E.multiply(inverseT);
        
        double lambda = new LUDecompositionImpl(EinverseT).getDeterminant();
        return lambda;
    }
    
    /**
     * Compute a Pillau Bartlett Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getPillauBartlettTrace(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Pillau-Bartlett Trace: hypothesis and error matrices must be square and same dimensions");
        
        RealMatrix T = H.add(E);
        RealMatrix inverseT = new LUDecompositionImpl(T).getSolver().getInverse();

        RealMatrix HinverseT = H.multiply(inverseT);
        
        return HinverseT.getTrace();
    }
    
    /**
     * Compute a Hotelling-Lawley Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getHotellingLawleyTrace(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Hotelling Lawley Trace: hypothesis and error matrices must be square and same dimensions");
        
        RealMatrix inverseE = new LUDecompositionImpl(E).getSolver().getInverse();
        RealMatrix HinverseE = H.multiply(inverseE);
        
        return HinverseE.getTrace();
    }
    
    /**
     * Compute a Univariate Approach to Repeated Measures statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getUnirep(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Unirep statistic: hypothesis and error matrices must be square and same dimensions");

        return H.getTrace() / E.getTrace();
    }
    
    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    private RealMatrix getHypothesisSumOfSquares(LinearModelPowerSampleSizeParameters params)
    {
        // convenience variables
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix B = params.getBeta();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix theta0 = params.getTheta();
        RealMatrix X = params.getDesign();
        
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(B.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(theta0);
        // get X'X invserse
        RealMatrix XtX = X.transpose().multiply(X);
        RealMatrix XtXInverse = new LUDecompositionImpl(XtX).getSolver().getInverse();
        
        // the middle term [C(X'X)-1C']-1
        RealMatrix cxxc = C.multiply(XtXInverse.multiply(C.transpose()));
        RealMatrix cxxcInverse = new LUDecompositionImpl(cxxc).getSolver().getInverse();
        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(cxxcInverse.multiply(thetaDiff));
        
        return hss;
        
    }
    
    /**
     * Calculate the sum of squares error matrix (the E matrix)
     * 
     * @param params matrices input by the user
     * @return sum o
     */
    private RealMatrix getErrorSumOfSquares(LinearModelPowerSampleSizeParameters params)
    {
        // get the rank of the design matrix, X
        int rankX = new SingularValueDecompositionImpl(params.getDesign()).getRank();
        // get the #rows in the design matrix
        int nX = params.getDesign().getRowDimension();
        
        RealMatrix U = params.getWithinSubjectContrast();
        return U.transpose().multiply(params.getSigmaError().multiply(U)).scalarMultiply(nX - rankX);
    }
}
