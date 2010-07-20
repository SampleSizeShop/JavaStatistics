package edu.cudenver.bios.power.glmm;

import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.MomentApproximationMethod;

public class GLMMTestHotellingLawley extends GLMMTest
{    
    public GLMMTestHotellingLawley(GLMMPowerParameters params)
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
        double a = (double) C.getRowDimension();
        // b = #columns in within subject contrast matrix
        double b = (double) U.getColumnDimension();
        // N = total number of subjects (rows in design matrix, X)
        double N = (double) X.getRowDimension();
        // r = rank of design matrix, X
        double r = (double) new SingularValueDecompositionImpl(X).getRank();
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * ((N - r) - b -1) + 2;
        }
        else
        {
            double t1 = (N - r) * (N - r) - (N - r) * (2 * b + 3) + b * (b + 3);
            double t2 = ((N - r) * (a  + b + 1) - (a + 2 * b + b * b - 1));
            df = 4 + (a * b + 2) * (t1/t2);
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
        RealMatrix B = params.getScaledBeta();
        
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
        if (type == DistributionType.DATA_ANALYSIS_NULL)
        {
            RealMatrix hypothesisSumOfSquares = getHypothesisSumOfSquares(params);
            RealMatrix errorSumOfSquares = getErrorSumOfSquares(params);
            double HLT = getHotellingLawleyTrace(hypothesisSumOfSquares, errorSumOfSquares);
            double ddf = getDenominatorDF(type);
            double ndf = getNumeratorDF(type);
            double b = params.getWithinSubjectContrast().getColumnDimension();
            RealMatrix X = params.getDesign();
            int r = new SingularValueDecompositionImpl(X).getRank();
            int N = X.getRowDimension();
            return HLT * (((N-r)-b-1)*ddf) / (ndf*(ddf-2));
        }
        else
        {
            return getNonCentrality(type) / getNumeratorDF(type);
        }
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
