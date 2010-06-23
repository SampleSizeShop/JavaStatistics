package edu.cudenver.bios.power.glmm;

import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;
import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters.MomentApproximationMethod;

public class GLMMTestHotellingLawley extends GLMMTest
{    
    public GLMMTestHotellingLawley(LinearModelPowerSampleSizeParameters params)
    {
        super(params);
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
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix X = params.getDesign();
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        
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
        
        double df = Double.NaN;
        if (type == DistributionType.DATA_ANALYSIS_NULL ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * ((N - r) - b -1) + 2;
        }
        else
        {
            df = (N - r) * (N - r) - (N - r) * (2 * b + 3) + b * (b + 3);
            df = df / ((N - r) * (a  + b + 1) - (a + 2 * b + b * b - 1));
            df = 4 + (a * b + 2) * df;
        }
        // TODO Auto-generated method stub
        return df;
    }

    @Override
    public double getNonCentrality(DistributionType type)
    {
        // calculate the hypothesis and error sum of squares matrices
        RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(params);
        RealMatrix errorSumOfSquares = getErrorSumOfSquares(params);
        
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix B = params.getBeta();
        
        // check if we are uni or multi variate
        double p = B.getColumnDimension();
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = U.getColumnDimension();
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
        
        if ((s == 1 && p > 1) ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                params.getMomentMethod() == MomentApproximationMethod.MCKEON_TWO_MOMENT_OMEGA_MULT)
        {
            RealMatrix X = params.getDesign();
            int r = new SingularValueDecompositionImpl(X).getRank();
            int N = X.getRowDimension();
            HLT *= ((double)(N - r)/(double)N);
            return N * s * HLT / s;
        }
        else
        {
            return getDenominatorDF(type) * HLT / s;
        }
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
    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        return a * b;
    }

    @Override
    public double getObservedF(DistributionType type)
    {                
        // a = #rows in between subject contrast matrix, C
        double a = params.getBetweenSubjectContrast().getRowDimension();
        // b = #columns in within subject contrast matrix, U
        double b = params.getWithinSubjectContrast().getColumnDimension();
                
        return getNonCentrality(type) / (a*b);
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
}
