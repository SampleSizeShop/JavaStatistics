package edu.cudenver.bios.power.glmm;

import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.MomentApproximationMethod;

public class GLMMTestPillaiBartlett extends GLMMTest
{
    public GLMMTestPillaiBartlett(GLMMPowerParameters params)
    {
        super(params);
    }   
    
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix X = params.getDesign();
        RealMatrix C = params.getBetweenSubjectContrast();
        RealMatrix U = params.getWithinSubjectContrast();
        
        // a = #rows in between subject contrast matrix, C
        double a = C.getRowDimension();
        // b = #columns in within subject contrast matrix
        double b = U.getColumnDimension();
        // N = total number of subjects (rows in design matrix, X)
        double N = X.getRowDimension();
        // r = rank of design matrix, X
        double r = new SingularValueDecompositionImpl(X).getRank();
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * ((N - r) - b + s);
        }
        else
        {
            double mu1= a * b / (N - r + a);
            double factor1 = (N - r + a - b) / (N - r + a - 1);
            double factor2 = (N - r) / (N - r + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((N - r + a)*(N - r + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * (m1 - m2) * (1 - m1) / denom;
        }

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
        
        double PB = getPillaiBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
        
        if ((s == 1 && p > 1) ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                params.getMomentMethod() == MomentApproximationMethod.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            RealMatrix X = params.getDesign();
            int N = X.getRowDimension();
            return N * s * PB / (s - PB);
        }
        else
        {
            return getDenominatorDF(type) * PB / (s - PB);
        }
    }

    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        double s = (a < b) ? a : b;  
        
        double df = Double.NaN;
        if (params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = a * b;
        }
        else
        {
            RealMatrix X = params.getDesign();
            // N = total number of subjects (rows in design matrix, X)
            double N = X.getRowDimension();
            // r = rank of design matrix, X
            double r = new SingularValueDecompositionImpl(X).getRank();
            double mu1= a * b / (N - r + a);
            double factor1 = (N - r + a - b) / (N - r + a - 1);
            double factor2 = (N - r) / (N - r + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((N - r + a)*(N - r + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * m1 * (m1 - m2) / denom;
        }
        
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
       // minimum of a and b dimensions
        double s = (a < b) ? a : b;  
        
        double association = 0.0;
        
        double PB = getPillaiBartlettTrace(hypothesisSumOfSquares, errorSumOfSquares);
        association = PB / s;
        double ddf = getDenominatorDF(type);
        double ndf = getNumeratorDF(type);

        return ((association) / ndf) / ((1 - association) / ddf);
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
        
        double a = params.getBetweenSubjectContrast().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        double s = (a < b) ? a : b;  
        double p = params.getBeta().getColumnDimension();
        
        RealMatrix adjustedH = H;
        if ((s == 1 && p > 1) ||
                params.getMomentMethod() == MomentApproximationMethod.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                params.getMomentMethod() == MomentApproximationMethod.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            RealMatrix X = params.getDesign();
            int r = new SingularValueDecompositionImpl(X).getRank();
            int N = X.getRowDimension();
            adjustedH = H.scalarMultiply(((double)(N - r)/(double)N));
        }
            
        RealMatrix T = adjustedH.add(E);
        RealMatrix inverseT = new LUDecompositionImpl(T).getSolver().getInverse();

        RealMatrix HinverseT = adjustedH.multiply(inverseT);
        
        return HinverseT.getTrace();
    }
}
