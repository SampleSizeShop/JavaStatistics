package edu.ucdenver.bios.statisticaltest;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.power.glmm.GLMMTest.DistributionType;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;

public class PillaiBartlettTraceTest extends StatisticalTest {

    private static final double TOLERANCE = 1.0E-15;
    // hypothesis sum of squares
    private RealMatrix hypothesisSumOfSquares;
    // error sum of squares
    private RealMatrix errorSumOfSquares;
    // value of the Pillai-Bartlett trace statistic
    private double PBT;

    // number of rows in the between participant contrast
    private double a;
    // number of columns in the within participant contrast
    private double b;
    // s = min(a,b)
    private double s;
    // total number of outcomes (i.e. columns of beta)
    private int p;
    // total sample size
    private int totalSampleSize;
    // nuE = designRank - total sample size
    private int nuE;
    
    // available f approximation methods
    private enum FApproximation
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

    // F approximation method
    FApproximation fMethod;
    
    /**
     * Constructor
     * @param hypothesisSumOfSquares
     * @param errorSumOfSquares
     * @param rowsBetweenContrast
     * @param columnsWithinContrast
     * @param totalOutcomes
     * @param designRank
     * @param totalSampleSize
     */
    public PillaiBartlettTraceTest(RealMatrix hypothesisSumOfSquares,
            RealMatrix errorSumOfSquares, 
            int rowsBetweenContrast, int columnsWithinContrast,
            int totalOutcomes, int designRank, int totalSampleSize) {

        if (!hypothesisSumOfSquares.isSquare() || 
                !errorSumOfSquares.isSquare() || 
                hypothesisSumOfSquares.getColumnDimension() != errorSumOfSquares.getRowDimension()) {
            throw new IllegalArgumentException("hypothesis and error matrices must be square and same dimensions");
        }

        this.hypothesisSumOfSquares = hypothesisSumOfSquares;
        this.errorSumOfSquares = errorSumOfSquares;
        this.a = rowsBetweenContrast;
        this.b = columnsWithinContrast;
        this.s = (a < b) ? a : b;  
        
        // calculate the value of the statistic
        this.PBT = getPillaiBartlettTrace();

        this.totalSampleSize = totalSampleSize;
        // calculate nuE
        this.nuE = totalSampleSize - designRank;
    }
    
    
    
    /**
     * Calculate the denominator degrees of freedom for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getDenominatorDF()
    {               
        double df = Double.NaN;
        if (fMethod == FApproximation.PILLAI_ONE_MOMENT ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * (nuE - b + s);
        }
        else
        {
            double mu1= a * b / (nuE + a);
            double factor1 = (nuE + a - b) / (nuE + a - 1);
            double factor2 = (nuE) / (nuE + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((nuE + a)*(nuE + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * (m1 - m2) * (1 - m1) / denom;
        }

        return df;
    }

    /**
     * Calculate the non-centrality parameter for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    private double getNonCentrality()
    {        
        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                fMethod == FApproximation.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            return totalSampleSize * s * PBT / (s - PBT);
        }
        else
        {
            return getDenominatorDF() * PBT / (s - PBT);
        }
    }

    /**
     * Calculate the numerator degrees of freedom for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getNumeratorDF()
    {
        double df = Double.NaN;
        if (
                fMethod == FApproximation.PILLAI_ONE_MOMENT ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = a * b;
        }
        else
        {
            double mu1= a * b / (nuE + a);
            double factor1 = (nuE + a - b) / (nuE + a - 1);
            double factor2 = (nuE) / (nuE + a + 2);
            double variance = 2 * a * b * factor1 * factor2 / ((nuE + a)*(nuE + a));
            double mu2 = variance + mu1 * mu1;
            double m1 = mu1 / s;
            double m2 = mu2 / (s*s);
            double denom = m2 - m1 * m1;
            df = 2 * m1 * (m1 - m2) / denom;
        }
        
        return df;
    }

    /**
     * Calculate the observed F for the PBT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return observed F
     * @throws IllegalArgumentException
     */
    private double getObservedF()
    {        
        double association = PBT / s;
        double ddf = getDenominatorDF();
        double ndf = getNumeratorDF();

        return ((association) / ndf) / ((1 - association) / ddf);
    }
    
    /**
     * Compute a Pillai Bartlett Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getPillaiBartlettTrace()
    {        
        RealMatrix adjustedH = hypothesisSumOfSquares;
        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                fMethod == FApproximation.MULLER_TWO_MOMENT_OMEGA_MULT)
        {
            adjustedH = hypothesisSumOfSquares.scalarMultiply(((double)(nuE)/(double)totalSampleSize));
        }
            
        RealMatrix T = adjustedH.add(errorSumOfSquares);
        RealMatrix inverseT = new LUDecomposition(T).getSolver().getInverse();

        RealMatrix HinverseT = adjustedH.multiply(inverseT);
        
        return HinverseT.getTrace();
    }
    
    
    
    
    @Override
    public AbstractRealDistribution getDataAnalysisNullDistribution() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public AbstractRealDistribution getPowerNullDistribution() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public AbstractRealDistribution getPowerAlternativeDistribution() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public double getStatistic() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double getPvalue() {
        // TODO Auto-generated method stub
        return 0;
    }

}
