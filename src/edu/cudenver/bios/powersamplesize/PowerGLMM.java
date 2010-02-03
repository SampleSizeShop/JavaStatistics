package edu.cudenver.bios.powersamplesize;

import java.io.FileWriter;

import jsc.distributions.FishersF;
import jsc.distributions.NoncentralFishersF;
import jsc.distributions.Normal;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;
import org.apache.commons.math.util.MathUtils;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.UnivariateCorrection;

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
    double unirepEpsilon = Double.NaN;
    double unirepEpsilonExpectedValue = Double.NaN;
    
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
        
        // get degrees of freedom for F distribution under null hypothesis
        // calculate numerator degrees of freedom
        double ndf = getNumeratorDF(powerParams, true);
        // calculate denominator df 
        double ddf = getDenominatorDF(powerParams, true);
        
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
            double fobs = getObservedF(simulatedParams, ddf);

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
        // unset the univariate epsilon - this is cached to avoid multiple calculations
        unirepEpsilon = Double.NaN;
        unirepEpsilonExpectedValue = Double.NaN;
        if (params.getTestStatistic() == TestStatistic.UNIREP)
            calculateUnirepCorrection(params);
        
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
        // get the degrees of freedom for the central F distribution under null hypothesis
        double nullNdf = getNumeratorDF(params, true);
        double nullDdf = getDenominatorDF(params, true);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params);
        
        // calculate the non-centrality parameter for the specified test statistic
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        double nonCentralityParam = a * b * getObservedF(params, nullDdf);
        // adjust for sphericity if using univariate approach to repeated measures test
        if (params.getTestStatistic() == TestStatistic.UNIREP) 
            nonCentralityParam *= unirepEpsilon;

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = getNumeratorDF(params, false);
        double altDdf = getDenominatorDF(params, false);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }
    
    private double getUnconditionalPower(LinearModelPowerSampleSizeParameters params)
    {
        // get the degrees of freedom for the central F distribution under null hypothesis
        double nullNdf = getNumeratorDF(params, true);
        double nullDdf = getDenominatorDF(params, true);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params);
        
        // determine unconditional power by integrating over all possible values 
        // of the non-centrality parameter
        // TODO: write me :-)
        
        return 0;
        
    }
    
    private double getQuantilePower(LinearModelPowerSampleSizeParameters params)
    {
        // get the degrees of freedom for the central F distribution under null hypothesis
        double nullNdf = getNumeratorDF(params, true);
        double nullDdf = getDenominatorDF(params, true);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params);
        
        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        PowerGLMMNonCentralityDistribution nonCentralityDist = 
            new PowerGLMMNonCentralityDistribution(params, false);
        double nonCentralityParam = nonCentralityDist.inverseCDF(params.getQuantile());
        
        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = getNumeratorDF(params, false);
        double altDdf = getDenominatorDF(params, false);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }
    
    /**
     * Calculate the critical F value for the specified test. 
     * 
     * @param ndf numerator degrees of freedom for central F
     * @param ddf denominator degrees of freedom for central F
     * @param params linear model inputs
     * 
     */
    private double getCriticalF(double ndf, double ddf, 
            LinearModelPowerSampleSizeParameters params)
    throws IllegalArgumentException
    {                                
        FishersF centralFDist = new FishersF(ndf, ddf);
        return centralFDist.inverseCdf(1 - params.getAlpha());
    }
    
    /**
     * Calculate the numerator degrees of freedom for the specified test, based on
     * whether the null or alternative hypothesis is assumed true.  Note, the null/alternative
     * degrees of freedom vary only for the corrected forms of the univariate approach to
     * repeated measures test.
     * 
     * @param powerParams input matrices
     * @param underNullHypothesis if true, returns degrees of freedom for central F distribution 
     * under the null hypothesis.  If false, returns degrees of freedom for non-central F under 
     * the alternative hypothesis.
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getNumeratorDF(LinearModelPowerSampleSizeParameters params, boolean underNullHypothesis)
    throws IllegalArgumentException
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        double df;
        
        if (params.getTestStatistic() == TestStatistic.UNIREP)
        {
            if (underNullHypothesis)
            {
                // unirep numerator df under null hypothesis
                switch(params.getUnivariateCorrection())
                {
                case GEISSER_GREENHOUSE:
                case HUYNH_FELDT:
                    // a * b * E[epsilon]
                    df = a * b * this.unirepEpsilonExpectedValue;
                    break;
                case BOX:
                    // a
                    df = a;
                    break;
                default: // uncorrected, a*b
                    df = a * b;
                    break;
                }
            }
            else
            {
                // unirep alternative hypotheses: a*b*epsilon 
                df = a * b * this.unirepEpsilon;
            }
        }
        else
        {
            // all other cases - same df used for both null and alternative
            // for wilks, HLT, PB
            df = a * b;
        }
        
        return df;
    }
    
    private void calculateUnirepCorrection(LinearModelPowerSampleSizeParameters params)
    {          
        // handle uncorrected unirep test
        if (params.getUnivariateCorrection() == UnivariateCorrection.NONE)
        {
            this.unirepEpsilon = 1;
            this.unirepEpsilonExpectedValue = 1;
        }
        else
        {
            // we are either not using the cached value or the correction epsilon
            // has not yet been calculated for this parameter set
            RealMatrix U = params.getWithinSubjectContrast();
            RealMatrix X = params.getDesign();
            int b = new SingularValueDecompositionImpl(U).getRank();
            int r = new SingularValueDecompositionImpl(U).getRank();
            RealMatrix sigmaStar = U.transpose().multiply(params.getSigmaError().multiply(U));
            double[] eigenValues = new EigenDecompositionImpl(sigmaStar, MathUtils.SAFE_MIN).getRealEigenvalues();
            double sumLambdaSquared = 0;
            double sumLambda = 0;
            for(double value: eigenValues)
            {
                sumLambda += value;
                sumLambdaSquared += value * value;
            }
            double epsilon = (sumLambda*sumLambda) / (b * (sumLambdaSquared));

            switch(params.getUnivariateCorrection())
            {
            case GEISSER_GREENHOUSE:
                unirepEpsilon = epsilon;
                break;
            case BOX:
                unirepEpsilon = epsilon;
                unirepEpsilonExpectedValue = 1;
                break;
            case HUYNH_FELDT:
                unirepEpsilon = 1;
                break;
            default:
                unirepEpsilon = 1;
            }
        }
    }
    
    /**
     * Calculate the denominator degrees of freedom for the specified test, based on
     * whether the null or alternative hypothesis is assumed true.  Note, the null/alternative
     * degrees of freedom vary only for the corrected forms of the univariate approach to
     * repeated measures test.
     * 
     * @param powerParams input matrices
     * @param underNullHypothesis if true, returns degrees of freedom for central F distribution 
     * under the null hypothesis.  If false, returns degrees of freedom for non-central F under 
     * the alternative hypothesis.
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getDenominatorDF(LinearModelPowerSampleSizeParameters powerParams,
            boolean underNullHypothesis)
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
        
        double df = -1;
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
        case PILLAI_BARTLETT_TRACE:
        {
            df = s * ((N - r) - b + s);
            break;
        }
        case UNIREP:
        {
            if (underNullHypothesis)
            {
                // unirep numerator df under null hypothesis
                switch(powerParams.getUnivariateCorrection())
                {
                case GEISSER_GREENHOUSE:
                case HUYNH_FELDT:
                    df = b*(N - r)*this.unirepEpsilonExpectedValue;
                    break;
                case BOX:
                    df = N - r;
                default: // uncorrected
                    df = b*(N - r);
                }
                df = b*(N - r)*this.unirepEpsilonExpectedValue;
            }
            else
            {
                // unirep under alternative hypothesis
                df = b*(N - r)*this.unirepEpsilon;
            }
            break;
        }
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
    private double getObservedF(LinearModelPowerSampleSizeParameters powerParams, double ddf)
    throws IllegalArgumentException
    {
        
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(powerParams);
        RealMatrix errorSumOfSquares = getErrorSumOfSquares(powerParams);
        
        RealMatrix C = powerParams.getBetweenSubjectContrast();
        RealMatrix U = powerParams.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
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
        case PILLAI_BARTLETT_TRACE:
            double PB = getPillaiBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
            association = PB / s;
            break;
        case UNIREP:
            double REP = getUnirep(hypothesisSumOfSquares, errorSumOfSquares);
            association = (REP / (1 + REP));
            break;
        case HOTELLING_LAWLEY_TRACE:
            double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
            association = (HLT/s) / (1 + (HLT/s));
            break;
        default:
            throw new IllegalArgumentException("Unknown statistic, cannot compute non-centrality parameter.");
        }
        
        double fobs = ((association) / (a*b)) / ((1 - association) / (double) ddf);
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
     * Compute a Pillai Bartlett Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getPillaiBartlettTrace(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Pillai-Bartlett Trace: hypothesis and error matrices must be square and same dimensions");
        
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
        int r = new SingularValueDecompositionImpl(params.getDesign()).getRank();
        // get the #rows in the design matrix
        int N = params.getDesign().getRowDimension();
        
        RealMatrix U = params.getWithinSubjectContrast();
        return U.transpose().multiply(params.getSigmaError().multiply(U)).scalarMultiply(N - r
                );
    }
}
