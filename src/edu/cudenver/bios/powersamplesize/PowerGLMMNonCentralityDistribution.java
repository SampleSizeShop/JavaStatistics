package edu.cudenver.bios.powersamplesize;

import jsc.distributions.FishersF;

import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.util.MathUtils;

import edu.cudenver.bios.matrix.ColumnMetaData;
import edu.cudenver.bios.matrix.EssenceMatrix;
import edu.cudenver.bios.matrix.ColumnMetaData.PredictorType;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.TestStatistic;

/**
 * Class representing the distribution of the non-centrality parameter in
 * the general linear multivariate model.  Used by the PowerGLMM class 
 * for computing unconditional and quantile power.
 * 
 * @see PowerGLMM
 * @author Sarah Kreidler
 */
public class PowerGLMMNonCentralityDistribution
{
    private static final double STARTING_NON_CENTRALITY = 100;    

    // intermediate forms 
    protected RealMatrix T1 = null;
    protected RealMatrix rootT1 = null;
    protected RealMatrix S = null;
    protected double H1;
    int qF;
    int a;
    int N;
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
    
    
    public PowerGLMMNonCentralityDistribution(LinearModelPowerSampleSizeParameters params, boolean exact)
    throws IllegalArgumentException
    {
        this.exact = exact;
        try
        {            
            // get design matrix for fixed parameters only
            RealMatrix F = params.getDesignEssence().getFullDesignFixed();
            qF = F.getColumnDimension();
            a = params.getBetweenSubjectContrast().getRowDimension();
            N = F.getRowDimension();
            // get fixed contrasts
            RealMatrix Cfixed = getFixedContrast(params);
            // build intermediate terms h1, S
            RealMatrix FtFinverse = 
                new LUDecompositionImpl(F.transpose().multiply(F)).getSolver().getInverse();
            RealMatrix P = Cfixed.multiply(FtFinverse).multiply(F.transpose());
            RealMatrix PPt = P.multiply(P.transpose());
            T1 = new LUDecompositionImpl(PPt).getSolver().getInverse();
            rootT1 = new CholeskyDecompositionImpl(T1).getL();
            RealMatrix theta0 = params.getTheta();
            RealMatrix C = params.getBetweenSubjectContrast();
            RealMatrix B = params.getBeta();
            RealMatrix U = params.getWithinSubjectContrast();
            // thetaHat = C * Beta * U
            RealMatrix thetaHat = C.multiply(B.multiply(U));
            // thetaHat - thetaNull.  Multiple by negative one to do subtraction
            RealMatrix thetaDiff = thetaHat.subtract(theta0);
            // TODO: specific to HLT or UNIREP
            RealMatrix sigmaStarInverse = getSigmaStarInverse(params);
            H1 = thetaDiff.transpose().multiply(T1).multiply(thetaDiff).multiply(sigmaStarInverse).getTrace();
            // matrix which represents the non-centrality parameter as a linear combination of chi-squared r.v.'s
            S = rootT1.transpose().multiply(thetaDiff).multiply(sigmaStarInverse).multiply(thetaDiff.transpose()).multiply(rootT1).scalarMultiply(1/H1);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e.getMessage());
        }

    }
    
    public double cdf(double w)
    {
        try
        {
            double b0 = 1 - w / H1;
            // we use the S matrix to generate the F-critical, numerical df's, and denominator df's 
            // for a central F distribution.  The resulting F distribution is used as an approximation
            // for the distribution of the non-centrality parameter
            // See formulas 18-21 and A8,A10 from Glueck & Muller (2003) for details
            double[] sEigenValues = new EigenDecompositionImpl(S, MathUtils.SAFE_MIN).getRealEigenvalues();
            // count the # of positive eigen values
            int sStar = 0;
            for(double value: sEigenValues) 
            {
                if (value > 0) sStar++;
            }
            double mz = 1; // TODO: calculate this correctly
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
            for(int k = 0; k < a; k++)
            {
                if (k == 0)
                {
                    // initial chi-square term is central (delta=0) with N-qf df, and lambda = b0
                    nu = N - qF;
                    lambda = b0;
                    delta = 0;
                }
                else if (k >= 1 && k <= sStar)
                {
                    // for k = 1 to sStar, chi-square term is non-central (delta = mz^2), 1 df,
                    // lambda = (b0 - kth eigen value of S)
                    nu = 1;
                    lambda = b0 - sEigenValues[k];
                    delta = mz * mz;
                }
                else
                {
                    // for k = sStar+1 to a, chi-sqaure term is non-central (delta = mz^2), 1 df,
                    // lambda = b0
                    nu = 1;
                    lambda = b0;
                    delta = mz * mz;
                }

                // accumulate terms
                if (lambda > 0)
                {
                    // positive terms
                    numPositive++;
                    m1Positive += lambda * (nu + delta);
                    m2Positive += lambda * lambda * 2* (nu + delta);
                }
                else if (lambda < 0)
                {
                    // negative terms - we take absolute value of lambda where needed
                    numNegative++;
                    m1Negative += -1 * lambda * (nu + delta);
                    m2Negative += lambda * lambda * 2* (nu + delta);
                }
                // TODO: what if lamba == 0?
            }

            // handle special cases

            // handle general case
            double nuStarPositive = 2 * (m1Positive * m1Positive) / m2Positive;
            double nuStarNegative = 2 * (m1Negative * m1Negative) / m2Negative;
            double lambdaStarPositive = m1Positive / (2 * m1Positive);
            double lambdaStarNegative =  m2Negative / (2 * m2Negative);
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
            double upperBound = getQuantileUpperBound(quantile);
            return solver.solve(quantFunc, 0, upperBound);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Failed to determine non-centrality quantile: " + e.getMessage());
        }
    }
    
    private double getQuantileUpperBound(double quantile)
    {
        double upperBound = STARTING_NON_CENTRALITY;
        
        for(double probNonCentral = 0.0; probNonCentral < quantile; upperBound *= 2)
        {
            probNonCentral = cdf(upperBound);
        }
        return upperBound;
    }
    
    private RealMatrix getFixedContrast(LinearModelPowerSampleSizeParameters params)
    {
        EssenceMatrix essenceX = params.getDesignEssence();
        RealMatrix C = params.getBetweenSubjectContrast();
        // determine which column contains the random predictor and include 
        // all columns representing fixed predictors
        // TODO: will need to modify this code for support of multiple random predictors
        int randCol = -1;
        int[] cols = new int[C.getColumnDimension() - essenceX.getRandomPredictorCount()];
        int colIdx = 0;
        for(int c = 0; c < C.getColumnDimension(); c++)
        {
            ColumnMetaData colMD = essenceX.getColumnMetaData(c);
            if (colMD.getPredictorType() == PredictorType.RANDOM)
            {
                randCol = c;
            }
            else if (colMD.getPredictorType() == PredictorType.FIXED)
            {
                cols[colIdx] = c;
                colIdx++;
            }
        }
        // count the number of rows which do not contain the random column
        int fixedRowCount = 0;
        for(int r = 0; r < C.getRowDimension(); r++)
        {
            if (C.getEntry(r, randCol) == 0) fixedRowCount++;
        }
        if (fixedRowCount == 0) throw new IllegalArgumentException("No fixed comparisons in contrast matrix");
        // now include the "fixed" rows of the between subject contrast in the final submatrix
        int[] rows = new int[fixedRowCount];
        int rowIdx = 0;
        for(int r = 0; r < C.getRowDimension(); r++)
        {
            if (C.getEntry(r, randCol) == 0) 
            {
                rows[rowIdx] = r;
                rowIdx++;
            }
        }

        return C.getSubMatrix(rows, cols); 
    }
    
    private RealMatrix getSigmaStarInverse(LinearModelPowerSampleSizeParameters params)
    {
        RealMatrix U = params.getWithinSubjectContrast();
        // sigma* = U'*sigmaE*U
        RealMatrix sigmaStar = U.transpose().multiply(params.getSigmaError()).multiply(U);
        
        if (params.getTestStatistic() == TestStatistic.UNIREP)
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
}
