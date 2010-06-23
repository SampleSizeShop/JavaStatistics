package edu.cudenver.bios.power.glmm;

import jsc.distributions.FishersF;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

/**
 * Class representing the distribution of the non-centrality parameter in
 * the general linear multivariate model.  Used by the PowerGLMM class 
 * for computing unconditional and quantile power.
 * 
 * @see PowerGLMM
 * @author Sarah Kreidler
 */
public class NonCentralityDistribution
{
    private static final double TOLERANCE = 0.000000000001;
    // intermediate forms 
    protected RealMatrix T1 = null;
    protected RealMatrix FT1 = null;
    protected RealMatrix S = null;
    protected RealMatrix mzSq = null;
    protected double H1;
    int qF;
    int a;
    int N;
    double[] sEigenValues;
    int sStar = 0;
    // indicates if an "exact" cdf should be calculated via Davie's algorithm or
    // with the Satterthwaite approximation from Glueck & Muller
    protected boolean exact;
    
    /**
     * Function calculating the difference between the probability of a target quantile 
     * and the  (used by the bisection solver from Apache Commons Math)
     * @see org.apache.commons.math.analysis.UnivariateRealFunction
     */
    private class NonCentralityQuantileFunction implements UnivariateRealFunction
    {
        protected double quantile;
        
        public NonCentralityQuantileFunction(double quantile)
        {
            this.quantile = quantile;
        }
        
        public double value(double n)
        {
            return cdf(n) - quantile;
        }
    }
    
    
    public NonCentralityDistribution(GLMMPowerParameters params, boolean exact)
    throws IllegalArgumentException
    {
        this.exact = exact;
        try
        {                    
            // TODO: need to calculate H0, need to adjust H1 for Unirep
            // get design matrix for fixed parameters only
            RealMatrix F = params.getDesignEssence().getFullDesignFixed();
            qF = F.getColumnDimension();
            a = params.getBetweenSubjectContrast().getRowDimension();
            N = F.getRowDimension();
            // get fixed contrasts
            RealMatrix Cfixed = getFixedContrast(params);
            RealMatrix CGaussian = getGaussianContrast(params);
            // build intermediate terms h1, S
            RealMatrix FtFinverse = 
                new LUDecompositionImpl(F.transpose().multiply(F)).getSolver().getInverse();
            RealMatrix P = Cfixed.multiply(FtFinverse).multiply(F.transpose());
            RealMatrix PPt = P.multiply(P.transpose());            
            T1 = new LUDecompositionImpl(PPt).getSolver().getInverse();
            FT1 = new CholeskyDecompositionImpl(T1).getL();
            // calculate theta difference
            RealMatrix theta0 = params.getTheta();
            RealMatrix C = params.getBetweenSubjectContrast();
            RealMatrix B = params.getScaledBeta();
            RealMatrix U = params.getWithinSubjectContrast();
            // thetaHat = C * Beta * U
            RealMatrix thetaHat = C.multiply(B.multiply(U));
            // thetaHat - thetaNull.  
            RealMatrix thetaDiff = thetaHat.subtract(theta0);
            // TODO: specific to HLT or UNIREP
            RealMatrix sigmaStarInverse = getSigmaStarInverse(params);
            RealMatrix H1matrix = thetaDiff.transpose().multiply(T1).multiply(thetaDiff).multiply(sigmaStarInverse);
            H1 = H1matrix.getTrace();
            // matrix which represents the non-centrality parameter as a linear combination of chi-squared r.v.'s
            S = FT1.transpose().multiply(thetaDiff).multiply(sigmaStarInverse).multiply(thetaDiff.transpose()).multiply(FT1).scalarMultiply(1/H1);
            // we use the S matrix to generate the F-critical, numerical df's, and denominator df's 
            // for a central F distribution.  The resulting F distribution is used as an approximation
            // for the distribution of the non-centrality parameter
            // See formulas 18-21 and A8,A10 from Glueck & Muller (2003) for details
            // TODO: move this into constructor since not dependent on w?
            EigenDecompositionImpl sEigenDecomp = new EigenDecompositionImpl(S, TOLERANCE);
            sEigenValues = sEigenDecomp.getRealEigenvalues();
            // count the # of positive eigen values
            for(double value: sEigenValues) 
            {
                if (value > 0) sStar++;
            }
            
            double stddevG = getStdDevG(params.getDesignEssence());
            mzSq = sEigenDecomp.getD().transpose().multiply(FT1.transpose()).multiply(CGaussian).scalarMultiply(1/stddevG);
            for(int i = 0; i < mzSq.getRowDimension(); i++)
            {
                for (int j = 0; j < mzSq.getColumnDimension(); j++)
                {
                    double entry = mzSq.getEntry(i, j);
                    mzSq.setEntry(i, j, entry*entry); // TODO: is there an apache function to do this?
                }
            }
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e.getMessage());
        }
    }
    
    public double cdf(double w)
    {
        if (H1 <= 0) return 0;
        if (H1 - w <= 0) return 1;
        try
        {
            double b0 = 1 - w / H1;
            double m1Positive = 0;
            double m1Negative = 0;
            double m2Positive = 0;
            double m2Negative = 0;
            //
            int numPositive = 0;
            int numNegative = 0;
            double nu;
            double delta;
            double lambda;
            
            // add in the first chi-squared term in the estimate of the non-centrality
            // (expressed as a sum of weighted chi-squared r.v.s)
            // initial chi-square term is central (delta=0) with N-qf df, and lambda = b0
            nu = N - qF;
            lambda = b0;
            delta = 0;
            // accumulate terms
            if (lambda > 0)
            {
                // positive terms
                numPositive++;
                m1Positive += lambda * (nu + delta);
                m2Positive += lambda * lambda * 2* (nu + 2*delta);
            }
            else if (lambda < 0)
            {
                // negative terms - we take absolute value of lambda where needed
                numNegative++;
                m1Negative += -1 * lambda * (nu + delta);
                m2Negative += lambda * lambda * 2* (nu + 2*delta);
            }
            
            // accumulate the remaining terms 
            for(int k = 0; k < sStar; k++)
            {
                if (k < sStar)
                {
                    // for k = 1 (well, 0 in java array terms and 1 in the paper) to sStar, chi-square term is 
                    // non-central (delta = mz^2), 1 df, lambda = (b0 - kth eigen value of S)
                    nu = 1;
                    lambda = b0 - sEigenValues[k];
                    delta = mzSq.getEntry(k, 0);
                }
                else
                {
                    // for k = sStar+1 to a, chi-sqaure term is non-central (delta = mz^2), 1 df,
                    // lambda = b0
                    nu = 1;
                    lambda = b0;
                    delta = mzSq.getEntry(k, 0);
                }
                // accumulate terms
                if (lambda > 0)
                {
                    // positive terms
                    numPositive++;
                    m1Positive += lambda * (nu + delta);
                    m2Positive += lambda * lambda * 2* (nu + 2*delta);
                }
                else if (lambda < 0)
                {
                    // negative terms - we take absolute value of lambda where needed
                    numNegative++;
                    m1Negative += -1 * lambda * (nu + delta);
                    m2Negative += lambda * lambda * 2* (nu + 2*delta);
                }
                // Note, we deliberately ignore terms for which lambda == 0
            }

            // handle special cases
            if (numNegative == 0) return 0;
            if (numPositive == 0) return 1;
            
            // handle general case
            double nuStarPositive = 2 * (m1Positive * m1Positive) / m2Positive;
            double nuStarNegative = 2 * (m1Negative * m1Negative) / m2Negative;
            double lambdaStarPositive = m2Positive / (2 * m1Positive);
            double lambdaStarNegative =  m2Negative / (2 * m1Negative);

            // create a central F to approximate the distribution of the non-centrality parameter
            FishersF centralFDist = new FishersF(nuStarPositive, nuStarNegative);
            // return power based on the non-central F
            return centralFDist.cdf((nuStarNegative*lambdaStarNegative)/(nuStarPositive*lambdaStarPositive));
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
    }
    
    public double inverseCDF(double quantile)
    {
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver solver = factory.newBisectionSolver();

        NonCentralityQuantileFunction quantFunc = new NonCentralityQuantileFunction(quantile);
        
        try
        {
            return solver.solve(quantFunc, 0, H1);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to determine non-centrality quantile: " + e.getMessage());
        }
    }
        
    /**
     * Returns the column representing the single random predictor if specified
     * (note, this function will require modification if we support multiple random predictors)
     * @param params
     * @return
     */
    private RealMatrix getGaussianContrast(GLMMPowerParameters params)
    {
        EssenceMatrix essenceX = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast();
        int randCol = -1;
        for(int c = 0; c < C.getColumnDimension(); c++)
        {
            ColumnMetaData colMD = essenceX.getColumnMetaData(c);
            if (colMD.getPredictorType() == PredictorType.RANDOM)
            {
                randCol = c;
                break; 
            }
        }
        if (randCol >= 0)
            return C.getColumnMatrix(randCol);
        else
            return null;      
    }
    
    private RealMatrix getFixedContrast(GLMMPowerParameters params)
    {
        EssenceMatrix essenceX = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast();
        // include all columns representing fixed predictors
        // TODO: will need to modify this code for support of multiple random predictors
        int[] cols = new int[C.getColumnDimension() - essenceX.getRandomPredictorCount()];
        int colIdx = 0;
        for(int c = 0; c < C.getColumnDimension(); c++)
        {
            ColumnMetaData colMD = essenceX.getColumnMetaData(c);
            if (colMD.getPredictorType() == PredictorType.FIXED)
            {
                cols[colIdx] = c;
                colIdx++;
            }
        }
        // include all rows
        int[] rows = new int[C.getRowDimension()];
        for(int r = 0; r < C.getRowDimension(); r++) rows[r] = r;

        return C.getSubMatrix(rows, cols); 
    }
    
    private RealMatrix getSigmaStarInverse(GLMMPowerParameters params)
    {
        RealMatrix U = params.getWithinSubjectContrast();
        // sigma* = U'*sigmaE*U
        RealMatrix sigmaStar = U.transpose().multiply(params.getScaledSigmaError()).multiply(U);
        
        if (params.getTest() == GLMMPowerParameters.Test.UNIREP)
        {
            int b = sigmaStar.getColumnDimension();
            // get discrepancy from sphericity for unirep test
            double sigmaStarTrace = sigmaStar.getTrace();
            double sigmaStarSquaredTrace = sigmaStar.multiply(sigmaStar).getTrace();
            double epsilon = (sigmaStarTrace*sigmaStarTrace) / ((double) b * sigmaStarSquaredTrace);
            RealMatrix identity = MatrixUtils.createRealIdentityMatrix(b);
            return identity.scalarMultiply((double) b * epsilon / sigmaStarTrace);
        }
        else 
        {
            // stat should only be HLT at this point  (exception is thrown by valdiateParams otherwise)
            return new LUDecompositionImpl(sigmaStar).getSolver().getInverse();
        }
    }
    
    private double getStdDevG(EssenceMatrix essence)
    throws IllegalArgumentException
    {
        for(int c = 0; c < essence.getColumnDimension(); c++)
        {
            ColumnMetaData colMD = essence.getColumnMetaData(c);
            if (colMD.getPredictorType() == PredictorType.RANDOM)
            {
                return Math.sqrt(colMD.getVariance());
            }
        }
        throw new IllegalArgumentException("Failed to calculate sigmaG - No random predictor");
    }

    public double getH1()
    {
        return H1;
    }   
    
}
