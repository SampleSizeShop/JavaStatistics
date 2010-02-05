package edu.cudenver.bios.powersamplesize;

import java.util.ArrayList;
import java.util.Arrays;

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
    
    // for the unirep test, the degrees of freedom change depending on two 
    // factors - data analysis (i.e. simulation or model fit) vs. power analysis,
    // and whether we need df for the F distribution under the null or the 
    // alternative hypothesis
    private enum DegreesOfFreedomType
    {
    	POWER_NULL,
    	POWER_ALTERNATIVE,
    	DATA_ANALYSIS_NULL
    };
    
    private class EigenValueMultiplicityPair
    {
    	public double eigenValue;
    	public double multiplicity;
    	
    	public EigenValueMultiplicityPair(double eigenValue, double multiplicity)
    	{
    		this.eigenValue = eigenValue;
    		this.multiplicity = multiplicity;
    	}
    };
    
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
        double ndf = getNumeratorDF(powerParams, DegreesOfFreedomType.DATA_ANALYSIS_NULL);
        // calculate denominator df 
        double ddf = getDenominatorDF(powerParams, DegreesOfFreedomType.DATA_ANALYSIS_NULL);
        
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
        double nullNdf = getNumeratorDF(params, DegreesOfFreedomType.POWER_NULL);
        double nullDdf = getDenominatorDF(params, DegreesOfFreedomType.POWER_NULL);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params.getAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        double nonCentralityParam = a * b * getObservedF(params, nullDdf);
        // adjust for sphericity if using univariate approach to repeated measures test
        if (params.getTestStatistic() == TestStatistic.UNIREP) 
            nonCentralityParam *= unirepEpsilon;

        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = getNumeratorDF(params, DegreesOfFreedomType.POWER_ALTERNATIVE);
        double altDdf = getDenominatorDF(params, DegreesOfFreedomType.POWER_ALTERNATIVE);
        // create a non-central F distribution (F distribution under the alternative hypothesis)
        NoncentralFishersF nonCentralFDist = new NoncentralFishersF(altNdf, altDdf, nonCentralityParam);
        // return power based on the non-central F
        return (1 - nonCentralFDist.cdf(Fcrit));  
    }
    
    private double getUnconditionalPower(LinearModelPowerSampleSizeParameters params)
    {
        // get the degrees of freedom for the central F distribution under null hypothesis
        double nullNdf = getNumeratorDF(params, DegreesOfFreedomType.POWER_NULL);
        double nullDdf = getDenominatorDF(params, DegreesOfFreedomType.POWER_NULL);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params.getAlpha());
        
        // determine unconditional power by integrating over all possible values 
        // of the non-centrality parameter
        // TODO: write me :-)
        
        return 0;
        
    }
    
    private double getQuantilePower(LinearModelPowerSampleSizeParameters params)
    {
        // get the degrees of freedom for the central F distribution under null hypothesis
        double nullNdf = getNumeratorDF(params, DegreesOfFreedomType.POWER_NULL);
        double nullDdf = getDenominatorDF(params, DegreesOfFreedomType.POWER_NULL);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = getCriticalF(nullNdf, nullDdf, params.getAlpha());
        
        // calculate the non-centrality parameter for the specified test statistic
        // For quantile power, we get the value from the distribution of the non-centrality
        // parameter which corresponds to the specified quantile
        PowerGLMMNonCentralityDistribution nonCentralityDist = 
            new PowerGLMMNonCentralityDistribution(params, false);
        double nonCentralityParam = nonCentralityDist.inverseCDF(params.getQuantile());
        
        // get the degrees of freedom for the non-central F under the alternative hypothesis
        // (these only change for the corrected Unirep test)
        double altNdf = getNumeratorDF(params, DegreesOfFreedomType.POWER_ALTERNATIVE);
        double altDdf = getDenominatorDF(params, DegreesOfFreedomType.POWER_ALTERNATIVE);
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
    private double getCriticalF(double ndf, double ddf, double alpha)
    throws IllegalArgumentException
    {                                
        FishersF centralFDist = new FishersF(ndf, ddf);
        return centralFDist.inverseCdf(1 - alpha);
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
    private double getNumeratorDF(LinearModelPowerSampleSizeParameters params, DegreesOfFreedomType dfType)
    throws IllegalArgumentException
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        // For Wilks, Hotelling-Lawley, and Pillai-Bartlett, the numerator df is always a* b
        double df = a * b;
        
        if (params.getTestStatistic() == TestStatistic.UNIREP)
        {
        	// for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        	// also if we are doing data analysis under the null hypothesis
    		switch(params.getUnivariateCorrection())
    		{
    		case GEISSER_GREENHOUSE:
    		case HUYNH_FELDT:
    		{
    			// GG, HF correction, we multiply the ndf by the epsilon estimate for
    			// power analysis (under alternative) and for data analysis.  For power under
    			// the null, we multiply by the expected value of the epsilon estimate
    			if (dfType == DegreesOfFreedomType.POWER_NULL)
    				df = a * b * this.unirepEpsilonExpectedValue;
    			else
    				df = a * b * this.unirepEpsilon;
    			break;
    		}
    		case BOX:
    			// for the conservative, or "Box" test, we adjust by epsilon only for
    			// power analysis under the alternative.  The ndf are the same for power
    			// under the null and for data analysis
    			if (dfType == DegreesOfFreedomType.POWER_ALTERNATIVE)
    				df = a * b * this.unirepEpsilon;
    			else
    				df = a;
    			break;
    		default: 
    			// in the uncorrected, we adjust by epsilon only for
    			// power analysis under the alternative.  The ndf are the same for power
    			// under the null and for data analysis
    			if (dfType == DegreesOfFreedomType.POWER_ALTERNATIVE)
    				df = a * b * this.unirepEpsilon;
    			else
    				df = a * b;
    			break;
    		}
        }
        
        return df;
    }
    
    private void calculateUnirepCorrection(LinearModelPowerSampleSizeParameters params)
    {          
    	RealMatrix U = params.getWithinSubjectContrast();
    	RealMatrix X = params.getDesign();
    	int b = new SingularValueDecompositionImpl(U).getRank();
    	int r = new SingularValueDecompositionImpl(U).getRank();
    	int N = X.getRowDimension();
    	// get the sigmaStar matrix: U' *sigmaError * U
    	RealMatrix sigmaStar = U.transpose().multiply(params.getSigmaError().multiply(U));
    	sigmaStar = sigmaStar.scalarMultiply(1/sigmaStar.getTrace()); // normalize sigma*

    	// get the eigen values of the normalized sigmaStar matrix
    	double[] eigenValues = new EigenDecompositionImpl(sigmaStar, MathUtils.SAFE_MIN).getRealEigenvalues();
    	if (eigenValues.length <= 0) throw new IllegalArgumentException("Failed to compute eigenvalues for sigma* matrix");
    	Arrays.sort(eigenValues); // put eigenvalues in ascending order

    	// calculate epsilon (correction for violation of sphericity)
    	// to avoid looping over the eigenvalues twice, we also calculate the multiplicity for distinct eigenvalues
    	
    	// list of distinct eigenvalues with multiplicity
    	ArrayList<EigenValueMultiplicityPair> distinctEigenValues = new ArrayList<EigenValueMultiplicityPair>();
    	// initialize values for the first eigen value
    	double first = eigenValues[0];
    	distinctEigenValues.add(new EigenValueMultiplicityPair(first, 1));
    	double sumLambda = first;
		double sumLambdaSquared = first * first;
		
    	// loop over remaining eigen values, saving distinct eigen values
    	for(int i = 1; i < eigenValues.length; i++)
    	{
    		double value = eigenValues[i];
    		// build the sum & sum of squares of eigen values
    		sumLambda += value;
    		sumLambdaSquared += value * value;
    		
    		// determine if this is a distinct eigen value and calculate multiplicity
			EigenValueMultiplicityPair prev = distinctEigenValues.get(distinctEigenValues.size()-1);
			if (Math.abs(prev.eigenValue - value) > MathUtils.SAFE_MIN)
			{
				// found new distinct eigen value
				distinctEigenValues.add(new EigenValueMultiplicityPair(value, 1));
			}
			else
			{
				// repeat of same eigenvalue, so  increment the multiplicity
				prev.multiplicity++;
			}
    	}
    	
    	// calculate estimate of epsilon (correction for violation of spehericity assumption)
    	unirepEpsilon = (sumLambda*sumLambda) / (b * (sumLambdaSquared));
    	// the Huynh-Feldt test uses a bias-corrected version of the epsilon estimate
    	if (params.getUnivariateCorrection() == UnivariateCorrection.HUYNH_FELDT)
    		unirepEpsilon = (N*b*unirepEpsilon - 2)/(b*(N - r - b*unirepEpsilon));

    	// For Geisser-Greenhouse and Huynh Feldt, we also need the expected value of the epsilon
    	// estimate.  Note, the estimates of epsilon represent different functions of the
    	// eigenvalues, so the resulting derivatives and expected values are specific to each test
    	if (params.getUnivariateCorrection() == UnivariateCorrection.GEISSER_GREENHOUSE)
    	{
    		// calculate the expected value of the epsilon estimate
    		// E[f(lambda)] = f(lambda) + g1 / (N - r)
    		// see Muller, Barton (1989) for details
    		double g1 = 0;
    		for(EigenValueMultiplicityPair evmI : distinctEigenValues)
    		{
    			double firstDerivative = 
    				((2 * sumLambda)/(b * sumLambdaSquared) - 
    					(2 * evmI.eigenValue * sumLambda * sumLambda)/(b*sumLambdaSquared*sumLambdaSquared));
    			
    			double secondDerivative = 
    				(2 / (b * sumLambdaSquared) - 
    					(8*evmI.eigenValue*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared) + 
    						(8*evmI.eigenValue*evmI.eigenValue*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared*sumLambdaSquared) - 
    							(2*sumLambda*sumLambda)/(b*sumLambdaSquared*sumLambdaSquared)); //TODO: finish
    			
    			// accumulate the first term of g1 (sum over distinct eigen vals of 1st derivative * eigen val ^2 * multiplicity)
    			g1 += secondDerivative * evmI.eigenValue * evmI.eigenValue * evmI.multiplicity;
    			// loop over elements not equal to current
    			for(EigenValueMultiplicityPair evmJ : distinctEigenValues)
    			{
    				if (evmI != evmJ)
    				{
    					// accumulate second term of g1
    					g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
    							(evmI.eigenValue - evmJ.eigenValue));
    				}
    			}
    		}
    		
    		this.unirepEpsilonExpectedValue = unirepEpsilon  + g1 / (N - r);
    	}
    	else if (params.getUnivariateCorrection() == UnivariateCorrection.HUYNH_FELDT)
    	{
    		// calculate the expected value of the epsilon estimate
    		// E[h(lambda)] = h(lambda) + g1 / (N - r)
    		// h(lambda) = h1(lambda) / (b*h2(lambda)
    		// see Muller, Barton (1989) for details
    		double h1 = N * sumLambda * sumLambda - 2 * sumLambdaSquared;
    		double h2 = (N - r) * sumLambdaSquared - (sumLambda * sumLambda);
    		double g1 = 0;
    		for(EigenValueMultiplicityPair evmI : distinctEigenValues)
    		{
    			// derivatives of sub-equations comprising epsilon estimator
    			double h1firstDerivative = (2 * N * sumLambda) - (4 *evmI.eigenValue); // TODO: finish
    			double h1secondDerivative = 2 * N - 4;
    				
    			double h2firstDerivative = (2 * (N - r) * evmI.eigenValue) - (2 * sumLambda); 
    			double h2SecondDerivative = 2 * (N - r) - 2;

    			// derivatives of estimate of epsilon
    			double firstDerivative = (h1firstDerivative / h2) - ((h1 * h2firstDerivative) / h2 * h2); 
    			double secondDerivative = 
    				((h1secondDerivative / h2) - (2 * h1firstDerivative * h2firstDerivative)/(h2 * h2) + 
    						(2 * h1 * h2firstDerivative * h2firstDerivative)/(h2 * h2 * h2) - (h1 * h2SecondDerivative)/(h2 * h2));
    			
    			// accumulate the first term of g1 (sum over distinct eigen vals of 1st derivative * eigen val ^2 * multiplicity)
    			g1 += secondDerivative * evmI.eigenValue * evmI.eigenValue * evmI.multiplicity;
    			// loop over elements not equal to current
    			for(EigenValueMultiplicityPair evmJ : distinctEigenValues)
    			{
    				if (evmI != evmJ)
    				{
    					// accumulate second term of g1
    					g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
    							(evmI.eigenValue - evmJ.eigenValue));
    				}
    			}
    		}
    		
    		this.unirepEpsilonExpectedValue = unirepEpsilon  + g1 / (N - r);
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
            DegreesOfFreedomType dfType)
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
        case HOTELLING_LAWLEY_TRACE:
        {
            df = s * ((N - r) - b -1) + 2;
            break;
        }
        case PILLAI_BARTLETT_TRACE:
        {
            df = s * ((N - r) - b + s);
            break;
        }
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
                df = (g*((N - r) - (b - a +1)/2)) - (a*b - 2)/2;
            }
            break;
        }
        case UNIREP:
        {
        	// for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        	// also if we are doing data analysis under the null hypothesis
        	switch(powerParams.getUnivariateCorrection())
        	{
        	case GEISSER_GREENHOUSE:
        	case HUYNH_FELDT:
        	{
        		// GG, HF correction, we multiply the ddf by the epsilon estimate for
        		// power analysis (under alternative) and for data analysis.  For power under
        		// the null, we multiply by the expected value of the epsilon estimate
        		if (dfType == DegreesOfFreedomType.POWER_NULL)
        			df = b*(N - r)*this.unirepEpsilonExpectedValue;
        		else
        			df = b*(N - r)*this.unirepEpsilon;
        		break;
        	}
        	case BOX:
        		// for the conservative, or "Box" test, we adjust by epsilon only for
        		// power analysis under the alternative.  The ddf are the same for power
        		// under the null and for data analysis
        		if (dfType == DegreesOfFreedomType.POWER_ALTERNATIVE)
        			df = b*(N - r) * this.unirepEpsilon;
        		else
        			df = (N - r);
        		break;
        	default: 
        		// in the uncorrected test, we adjust by epsilon only for
        		// power analysis under the alternative.  The ddf are the same for power
        		// under the null and for data analysis
        		if (dfType == DegreesOfFreedomType.POWER_ALTERNATIVE)
        			df = b*(N - r) * this.unirepEpsilon;
        		else
                    df = b*(N - r);
        		break;
        	}
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
        
        return ((association) / (a*b)) / ((1 - association) / ddf);
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
