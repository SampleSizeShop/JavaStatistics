package edu.cudenver.bios.power.glmm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;
import org.apache.commons.math.util.MathUtils;

import edu.cudenver.bios.matrix.GramSchmidtOrthonormalization;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMTestUnivariateRepeatedMeasures extends GLMMTest
{    
    protected static final double TOLERANCE = 0.000000000001;
    protected double unirepEpsilon = Double.NaN;
    
    protected class EigenValueMultiplicityPair
    {
        public double eigenValue;
        public double multiplicity;
        
        public EigenValueMultiplicityPair(double eigenValue, double multiplicity)
        {
            this.eigenValue = eigenValue;
            this.multiplicity = multiplicity;
        }
    };
    
    public GLMMTestUnivariateRepeatedMeasures(GLMMPowerParameters params)
    {
        super(params);
        
        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();
        
        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection();
    }
    
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix X = params.getDesign();
        RealMatrix U = params.getWithinSubjectContrast();
        
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();
        // N = total number of subjects (rows in design matrix, X)
        int N = X.getRowDimension();
        // r = rank of design matrix, X
        int r = new SingularValueDecompositionImpl(X).getRank();
        
        double df = Double.NaN;
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected test, we adjust by epsilon only for
        // power analysis under the alternative.  The ddf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = b*(N - r) * this.unirepEpsilon;
        else
            df = b*(N - r);
        
        return df;
    }

    @Override
    public double getNonCentrality(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*unirepEpsilon;
    }

    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        double df = Double.NaN;
        // for the unirep test, the degrees of freedom change for power under the null vs alternative, and
        // also if we are doing data analysis under the null hypothesis

        // in the uncorrected, we adjust by epsilon only for
        // power analysis under the alternative.  The ndf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = a * b * this.unirepEpsilon;
        else
            df = a * b;

        return df;
    }

    @Override
    public double getObservedF(DistributionType type)
    {
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(params);
        RealMatrix errorSumOfSquares = getErrorSumOfSquares(params);
        
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
        
        double association = 0.0;
        
        double REP = getUnirep(hypothesisSumOfSquares, errorSumOfSquares);
        association = (REP / (1 + REP));
        
        double ddf = getDenominatorDF(type);
        double ndf = getNumeratorDF(type);
        return ((association) / ndf) / ((1 - association) / ddf);
    }
    
    /**
     * Compute a Univariate Approach to Repeated Measures statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    protected double getUnirep(RealMatrix H, RealMatrix E)
    {
        if (!H.isSquare() || !E.isSquare() || H.getColumnDimension() != E.getRowDimension())
            throw new InvalidMatrixException("Failed to compute Unirep statistic: hypothesis and error matrices must be square and same dimensions");

        return H.getTrace() / E.getTrace();
    }
    
    private void calculateUnirepCorrection()
    {          
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix X = params.getDesign();
        int b = new SingularValueDecompositionImpl(U).getRank();
        int r = new SingularValueDecompositionImpl(X).getRank();
        int N = X.getRowDimension();
        // get the sigmaStar matrix: U' *sigmaError * U
        RealMatrix sigmaStar = U.transpose().multiply(params.getScaledSigmaError().multiply(U));
        // ensure symmetry
        sigmaStar = sigmaStar.add(sigmaStar.transpose()).scalarMultiply(0.5); 
        // normalize
        sigmaStar = sigmaStar.scalarMultiply(1/sigmaStar.getTrace());
        
        // get the eigen values of the normalized sigmaStar matrix
        double[] eigenValues = new EigenDecompositionImpl(sigmaStar, TOLERANCE).getRealEigenvalues();
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
            if (Math.abs(prev.eigenValue - value) > TOLERANCE)
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
    }
    
    protected void createOrthonormalU()
    {
        RealMatrix U = params.getWithinSubjectContrast();
        
        RealMatrix UtU = U.transpose().multiply(U);
        double upperLeft = UtU.getEntry(0, 0);
        if (upperLeft != 0) UtU = UtU.scalarMultiply(1/upperLeft);
        
        RealMatrix diffFromIdentity = UtU.subtract(MatrixUtils.createRealIdentityMatrix(UtU.getRowDimension()));
        
        // get maximum value in U'U
        double maxValue = Double.NEGATIVE_INFINITY;
        for(int r = 0; r < diffFromIdentity.getRowDimension(); r++)
        {
            for(int c = 0; c < diffFromIdentity.getColumnDimension(); c++)
            {
                double entryVal = diffFromIdentity.getEntry(r, c);
                if (entryVal > maxValue) maxValue = entryVal;
            }
        }
        
        if (maxValue > MathUtils.SAFE_MIN)
        {
            // U matrix deviates from Identity, so create one that is orthonormal
            RealMatrix Utmp = new GramSchmidtOrthonormalization(U).getQ();
            params.setWithinSubjectContrast(Utmp);
        }
    }
}
