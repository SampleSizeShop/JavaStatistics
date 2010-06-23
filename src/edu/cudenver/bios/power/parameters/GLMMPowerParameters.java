package edu.cudenver.bios.power.parameters;

import java.util.ArrayList;
import java.util.ListIterator;

import org.apache.commons.math.linear.RealMatrix;
import edu.cudenver.bios.matrix.EssenceMatrix;

/**
 * Container class for matrix inputs for general linear model power calculations.
 * <ul>
 * <li>Design matrix - essence matrix for the model design
 * <li>Beta matrix - estimated regression coefficients matrix
 * <li>Sigma matrix - estimated errors matrix
 * <li>Theta0 matrix - estimated null hypothesis matrix 
 * <li>Between subject contrast matrix - defines comparisons between subjects
 * <li>Within subject contrast matrix - defines comparisons within subjects
 * (may be left null for univariate special case)
 * </ul>
 * @author kreidles
 *
 */
public class GLMMPowerParameters extends PowerParameters
{
	// the type of statistical test to use
	public enum Test 
	{
	    NONE,
		UNIREP,
		UNIREP_BOX,
		UNIREP_GEISSER_GREENHOUSE,
		UNIREP_HUYNH_FELDT,
		WILKS_LAMBDA,
		PILLAI_BARTLETT_TRACE,
		HOTELLING_LAWLEY_TRACE
	};

	public enum PowerMethod
	{
		CONDITIONAL_POWER,
		UNCONDITIONAL_POWER,
		QUANTILE_POWER
	};

	public enum UnivariateCdf
	{
		MULLER_BARTON_APPROX,
		MULLER_EDWARDS_TAYLOR_APPROX,
		MULLER_EDWARDS_TAYLOR_EXACT,
		MULLER_EDWARDS_TAYLOR_EXACT_APPROX
	};

	public enum MomentApproximationMethod
	{
		NONE,
		PILLAI_ONE_MOMENT,
		PILLAI_ONE_MOMENT_OMEGA_MULT,
		MCKEON_TWO_MOMENT,
		MCKEON_TWO_MOMENT_OMEGA_MULT,
		MULLER_TWO_MOMENT,
		MULLER_TWO_MOMENT_OMEGA_MULT,
		RAO_TWO_MOMENT,
		RAO_TWO_MOMENT_OMEGA_MULT
	};

	Test test = Test.NONE;

	RealMatrix beta = null;
	RealMatrix betaScaled = null;
	PeekableList<Double> betaScaleList = new PeekableList<Double>();
	
	// used if only fixed predictors
	RealMatrix sigmaError = null;
	RealMatrix sigmaErrorScaled = null;
	PeekableList<Double> sigmaScaleList = new PeekableList<Double>();
	
	/* variance/ covariances required for random gaussian covariates */
	RealMatrix sigmaOutcomeGaussianRandom = null;
	RealMatrix sigmaGaussianRandom = null;
	RealMatrix sigmaOutcome = null;

	RealMatrix betweenSubjectContrast = null;

	RealMatrix theta = null;

	RealMatrix withinSubjectContrast = null;

	RealMatrix design = null;

	EssenceMatrix designEssence = null;

	PowerMethod powerMethod = PowerMethod.CONDITIONAL_POWER;
	double quantile = 0.50;

	UnivariateCdf univariateCdf = UnivariateCdf.MULLER_EDWARDS_TAYLOR_APPROX;

	MomentApproximationMethod momentMethod =
		MomentApproximationMethod.NONE;

	/**
	 * Constructor.  Creates an empty set of linear model power parameters
	 */
	public GLMMPowerParameters() {}

	public Test getTest()
	{
		return test;
	}

	public void setTest(Test test)
	{
	    switch (test)
	    {
	    case HOTELLING_LAWLEY_TRACE:
	        momentMethod = MomentApproximationMethod.MCKEON_TWO_MOMENT_OMEGA_MULT;
	        break;
	    case PILLAI_BARTLETT_TRACE:
	        momentMethod = MomentApproximationMethod.MULLER_TWO_MOMENT;
	        break;
	    case WILKS_LAMBDA:
	        momentMethod = MomentApproximationMethod.RAO_TWO_MOMENT_OMEGA_MULT;
	        break;
	    default:
	        momentMethod = MomentApproximationMethod.NONE;
	    }

		this.test = test;
	}

	public RealMatrix getBeta()
	{
		return beta;
	}

	public void setBeta(RealMatrix beta)
	{
		this.beta = beta;
	}

	public RealMatrix getSigmaError()
	{
		return sigmaError;
	}

	public void setSigmaError(RealMatrix sigmaError)
	{
		this.sigmaError = sigmaError;
	}

	public RealMatrix getSigmaOutcomeGaussianRandom()
	{
		return sigmaOutcomeGaussianRandom;
	}

	public void setSigmaOutcomeGaussianRandom(RealMatrix sigmaOutcomeGaussianRandom)
	{
		this.sigmaOutcomeGaussianRandom = sigmaOutcomeGaussianRandom;
	}

	public RealMatrix getSigmaGaussianRandom()
	{
		return sigmaGaussianRandom;
	}

	public void setSigmaGaussianRandom(RealMatrix sigmaGaussianRandom)
	{
		this.sigmaGaussianRandom = sigmaGaussianRandom;
	}

	public RealMatrix getSigmaOutcome()
	{
		return sigmaOutcome;
	}

	public void setSigmaOutcome(RealMatrix sigmaOutcome)
	{
		this.sigmaOutcome = sigmaOutcome;
	}

	public RealMatrix getBetweenSubjectContrast()
	{
		return betweenSubjectContrast;
	}

	public void setBetweenSubjectContrast(RealMatrix betweenSubjectContrast)
	{
		this.betweenSubjectContrast = betweenSubjectContrast;
	}

	public RealMatrix getTheta()
	{
		return theta;
	}

	public void setTheta(RealMatrix theta)
	{
		this.theta = theta;
	}

	public RealMatrix getWithinSubjectContrast()
	{
		return withinSubjectContrast;
	}

	public void setWithinSubjectContrast(RealMatrix withinSubjectContrast)
	{
		this.withinSubjectContrast = withinSubjectContrast;
	}

	public RealMatrix getDesign()
	{
		if (design == null && designEssence != null)
		{
			// if only a design matrix is specified, retrieve the full design matrix
			this.design = designEssence.getFullDesignMatrix();
			// TODO: what to do if this is too big to cache?
		}
		return design;
	}

	public void setDesign(RealMatrix design)
	{
		this.design = design;
	}

	public EssenceMatrix getDesignEssence()
	{
		return designEssence;
	}

	public void setDesignEssence(EssenceMatrix designEssence)
	{
		this.designEssence = designEssence;
	}

	public PowerMethod getPowerMethod()
	{
		return powerMethod;
	}

	public void setPowerMethod(
			PowerMethod powerMethod)
	{
		this.powerMethod = powerMethod;
	}

	public double getQuantile()
	{
		return quantile;
	}

	public void setQuantile(double quantile)
	{
		this.quantile = quantile;
	}

	public MomentApproximationMethod getMomentMethod()
	{
		return momentMethod;
	}

	public void setMomentMethod(MomentApproximationMethod momentMethod)
	{
		this.momentMethod = momentMethod;
	}

	public UnivariateCdf getUnivariateCdf()
	{
		return univariateCdf;
	}

	public void setUnivariateCdf(UnivariateCdf univariateCdf)
	{
		this.univariateCdf = univariateCdf;
	}

    /**
     * Add a scale factor for the sigma matrix
     */
    public void addSigmaScale(double scale)
    {
    	sigmaScaleList.add(new Double(scale));    	
    }
    
    /**
     * Add a scale factor for the beta matrix
     */
    public void addBetaScale(double scale)
    {
    	betaScaleList.add(new Double(scale));    	
    }    
	
    public Double getFirstBetaScale()
    {
        Double betaScale = betaScaleList.first();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue());
        return betaScale;
    }
    
    public Double getNextBetaScale()
    {
        Double betaScale = betaScaleList.next();
        if (betaScale != null)
            betaScaled = beta.scalarMultiply(betaScale.doubleValue());
        return betaScale;
    }
    
    public Double getCurrentBetaScale()
    {
        return betaScaleList.current();
    }
	
	public RealMatrix getScaledBeta()
	{
	    return (betaScaled != null ? betaScaled : beta);
	}
	
	
    public Double getFirstSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.first();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;        
    }
    
    public Double getNextSigmaScale()
    {
        Double sigmaScale = sigmaScaleList.next();
        if (sigmaScale != null)
            sigmaErrorScaled = sigmaError.scalarMultiply(sigmaScale.doubleValue());
        return sigmaScale;    
    }
    
    public Double getCurrentSigmaScale()
    {
        return sigmaScaleList.current();
    }
    
    public RealMatrix getScaledSigmaError()
    {
        return (sigmaErrorScaled != null ? sigmaErrorScaled : sigmaError);
    }

}


