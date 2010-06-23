package edu.cudenver.bios.power;

import java.util.ArrayList;
import java.util.List;

import jsc.distributions.NoncentralFishersF;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.Test;
import edu.cudenver.bios.power.glmm.GLMMTest;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.NonCentralityDistribution;

public class GLMMPowerCalculator implements PowerCalculator
{
	
    private static final int STARTING_SAMPLE_SIZE = 1000;
    private static final int MIN_SAMPLE_SIZE =  2; // need df > 0
    
    /**
     * function to be used with apache's built-in bisection solver
     */
    private class SampleSizeFunction implements UnivariateRealFunction
    {
    	GLMMPowerParameters params;
        
        public SampleSizeFunction(GLMMPowerParameters params)
        {
            this.params = params;
        }
        
        public double value(double n)
        {
            try
            {
//            	GLMMPowerParameters tmpParams = 
//                    new GLMMPowerParameters(params);
//                RealMatrix design = tmpParams.getDesignEssence().getFullDesignMatrix((int) n);
//                tmpParams.setDesign(design);
//                tmpParams.setSampleSize((int) n);
//                double calculatedPower = powerGLMM.getCalculatedPower(tmpParams);
//                return tmpParams.getPower() - calculatedPower;
                return 0;
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
     * Class passed into Apache's TrapezoidIntegrator function to compute
     * unconditional power
     */
    private class UnconditionalPowerIntegrand implements UnivariateRealFunction
    {
        protected NonCentralityDistribution nonCentralityDist;
        protected double Fcrit;
        protected double ndf;
        protected double ddf;

        public UnconditionalPowerIntegrand(NonCentralityDistribution nonCentralityDist,
                double Fcrit, double ndf, double ddf)
        {
            this.nonCentralityDist = nonCentralityDist;
            this.Fcrit = Fcrit;
            this.ndf = ndf;
            this.ddf = ddf;
        }
        
        public double value(double t)
        {
            NoncentralFishersF FdistTerm1 = new NoncentralFishersF(ndf, ddf, t);
            NoncentralFishersF FdistTerm2 = new NoncentralFishersF(ndf+2, ddf, t);

            return nonCentralityDist.cdf(t)*(FdistTerm1.cdf(Fcrit) - FdistTerm2.cdf((Fcrit*ndf)/(ndf+2)));
        }
    }
	
	/********* public methods for the power API ************/
	
	/**
	 * 
	 */
	@Override
	public List<Power> getPower(PowerParameters powerParams)
	{
        GLMMPowerParameters params = (GLMMPowerParameters) powerParams;
        // make sure all of the matrix inputs are appropriate
        validateMatrices(params);
        
        // precalculate any computationally expensive matrices/constants, 
        // update the parameters as needed - used for random covariates
        initialize(params);
        
        // list of power results
        ArrayList<Power> results = new ArrayList<Power>();
        
        // calculate the power for either one or two tails
        for(Double alpha = params.getFirstAlpha(); alpha != null;
        alpha = params.getNextAlpha())
        {
            for(Double betaScale = params.getFirstBetaScale(); betaScale != null;
            betaScale = params.getNextBetaScale())
            {
                for(Double sigmaScale = params.getFirstSigmaScale(); sigmaScale != null;
                sigmaScale = params.getNextSigmaScale())
                {
                    for(Integer sampleSize = params.getFirstSampleSize(); sampleSize != null; 
                    sampleSize = params.getNextSampleSize())
                    {
        			    // calculate the power
        			    double power = Double.NaN;
        				switch (params.getPowerMethod())
        		        {
        		        case QUANTILE_POWER:
        		            power = getQuantilePower(params);
        		        case UNCONDITIONAL_POWER:
        		            power = getUnconditionalPower(params);
        		        case CONDITIONAL_POWER:
        		        default:
        		            power = getConditionalPower(params);
        		        }
        				
        				// store the power result
                        results.add(new GLMMPower(alpha, power, power, sampleSize, betaScale, sigmaScale));
        			}
        		}
        	}
        }
        
        return results;
	}

	@Override
	public List<Power> getSampleSize(PowerParameters params)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Power> getSimulatedPower(PowerParameters params, int iterations)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Power> getDetectableDifference(PowerParameters params)
	{
		// TODO Auto-generated method stub
		return null;
	}
	
    private void initialize(GLMMPowerParameters params)
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
	
	protected void validateMatrices(GLMMPowerParameters params) throws IllegalArgumentException
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
//        if (params.getAlpha() <= 0 || params.getAlpha() >= 1)
//            throw new IllegalArgumentException("Alpha must be between 0 and 1.  Invalid value: " + params.getAlpha());
//        // make sure sample size > 1
//        if (params.getSampleSize() <= 1)
//            throw new IllegalArgumentException("Sample size must be greater than 1.  Invalid value: " + params.getSampleSize());
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
            if (params.getTest() != Test.HOTELLING_LAWLEY_TRACE &&
                    params.getTest() != Test.UNIREP)
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
	
    private double getConditionalPower(GLMMPowerParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
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
    
    private double getUnconditionalPower(GLMMPowerParameters params)
    throws IllegalArgumentException
    {  		
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
        // get the distribution of the noncentrality parameter
        NonCentralityDistribution nonCentralityDist = new NonCentralityDistribution(params, false);
        double ndf = glmmTest.getNumeratorDF(GLMMTest.DistributionType.POWER_NULL);
        double ddf = glmmTest.getDenominatorDF(GLMMTest.DistributionType.POWER_NULL);
        double h1 = nonCentralityDist.getH1();
        // integrate over all value of non-centrality parameter from h0 to h1
        UnconditionalPowerIntegrand integrand = 
            new UnconditionalPowerIntegrand(nonCentralityDist, Fcrit, ndf, ddf);
        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        try
        {
            // create a noncentral F dist with non-centrality of H1
            NoncentralFishersF fdist = new NoncentralFishersF(ndf, ddf, h1);
            double integralResult = integrator.integrate(integrand, 0, h1);
            
            return (1 - fdist.cdf(Fcrit) - 0.5*integralResult);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to integrate over non-centrality parameter: " + e.getMessage());
        }        
    }
    
    private double getQuantilePower(GLMMPowerParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, params.getCurrentAlpha());
        
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
