package edu.cudenver.bios.power.glmm;

import jsc.distributions.FishersF;

import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public abstract class GLMMTest
{    
    // for the unirep test, the degrees of freedom change depending on two 
    // factors - data analysis (i.e. simulation or model fit) vs. power analysis,
    // and whether we need df for the F distribution under the null or the 
    // alternative hypothesis
    public enum DistributionType
    {
        POWER_NULL,
        POWER_ALTERNATIVE,
        DATA_ANALYSIS_NULL
    };
    
    // store incoming parameters
    protected GLMMPowerParameters params;

    // cache some common operations like rank, sample size, etc.
    protected double N; // total sample size (row dimension of design matrix)
    protected double r; // rank of the design matrix
    
    public GLMMTest(GLMMPowerParameters params)
    {
        this.params = params;
        RealMatrix X = params.getDesign();
        N = X.getRowDimension();
        r = params.getDesignRank();
    }

    /**
     * Calculate the critical F value for the specified test. 
     * 
     * @param ndf numerator degrees of freedom for central F
     * @param ddf denominator degrees of freedom for central F
     * @param params linear model inputs
     * 
     */
    public double getCriticalF(DistributionType type, double alpha)
    throws IllegalArgumentException
    {                               
        
        double ndf = getNumeratorDF(type);
        double ddf = getDenominatorDF(type);

        FishersF centralFDist = new FishersF(ndf, ddf);
        return centralFDist.inverseCdf(1 - alpha);
    }
    
    abstract public double getNumeratorDF(DistributionType type);
    
    abstract public double getDenominatorDF(DistributionType type);
    
    abstract public double getObservedF(DistributionType type);
    
    abstract public double getNonCentrality(DistributionType type);
    
    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    protected RealMatrix getHypothesisSumOfSquares(GLMMPowerParameters params)
    {
        // convenience variables
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix B = params.getScaledBeta();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix theta0 = params.getTheta();
        
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(B.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(theta0);
        // get X'X invserse
        RealMatrix XtXInverse = params.getXtXInverse();
        
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
    protected RealMatrix getErrorSumOfSquares(GLMMPowerParameters params)
    {        
        RealMatrix U = params.getWithinSubjectContrast();
        return U.transpose().multiply(params.getScaledSigmaError().multiply(U)).scalarMultiply(N-r);
    }    
    
}
