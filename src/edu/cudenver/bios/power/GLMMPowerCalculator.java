package edu.cudenver.bios.power;

import java.util.ArrayList;
import java.util.List;

import jsc.distributions.NoncentralFishersF;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.OneSampleStudentsTPowerParameters;
import edu.cudenver.bios.power.parameters.PowerParameters;
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
        for(Double alpha: params.getAlphaList())
        {
        	for(Double betaScale: params.getBetaScaleList())
        	{
        		for(Double sigmaScale: params.getSigmaScaleList())
        		{
        			for(Integer sampleSize: params.getSampleSizeList())
        			{
        				switch (params.getPowerMethod())
        		        {
        		        case QUANTILE_POWER:
        		        {
        		        	results.add(getQuantilePower(alpha.doubleValue(), betaScale.doubleValue(), 
        		        			sigmaScale.doubleValue(), sampleSize.intValue(), params));
        		        }
        		        case UNCONDITIONAL_POWER:
        		        {
        		        	results.add(getUnconditionalPower(alpha.doubleValue(), betaScale.doubleValue(), 
        		        			sigmaScale.doubleValue(), sampleSize.intValue(), params));
        		        }
        		        case CONDITIONAL_POWER:
        		        default:
        		        {
        		        	results.add(getConditionalPower(alpha.doubleValue(), betaScale.doubleValue(), 
        		        			sigmaScale.doubleValue(), sampleSize.intValue(), params));
        		        }
        		        }
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
		
	}
	
    private GLMMPower getConditionalPower(double alpha, double betaScale, double sigmaScale, 
    		int sampleSize, GLMMPowerParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
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
        double power = 1 - nonCentralFDist.cdf(Fcrit);
        return new GLMMPower(alpha, power, power, sampleSize, betaScale, sigmaScale);
    }
    
    private GLMMPower getUnconditionalPower(double alpha, double betaScale, double sigmaScale, 
    		int sampleSize, GLMMPowerParameters params)
    throws IllegalArgumentException
    {  		
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);

        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
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
            
            double power = 1 - fdist.cdf(Fcrit) - 0.5*integralResult;
            return new GLMMPower(alpha, power, power, sampleSize, betaScale, sigmaScale);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to integrate over non-centrality parameter: " + e.getMessage());
        }        
    }
    
    private GLMMPower getQuantilePower(double alpha, double betaScale, double sigmaScale, 
    		int sampleSize, GLMMPowerParameters params)
    {
        GLMMTest glmmTest = GLMMTestFactory.createGLMMTest(params);
        // get the approximate critical F value (central F) under the null hypothesis
        double Fcrit = glmmTest.getCriticalF(GLMMTest.DistributionType.POWER_NULL, alpha);
        
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
        double power = 1 - nonCentralFDist.cdf(Fcrit);  
        
        return new GLMMPower(alpha, power, power, sampleSize, betaScale, sigmaScale);
    }

}
