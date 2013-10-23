package edu.ucdenver.bios.statisticaltest;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.power.glmm.GLMMTest.DistributionType;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;

public class WilksLambdaTest extends StatisticalTest {
    
    private static final double TOLERANCE = 1.0E-15;
    // hypothesis sum of squares
    private RealMatrix hypothesisSumOfSquares;
    // error sum of squares
    private RealMatrix errorSumOfSquares;
    // value of the HLT statistic
    private double lambda;

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
    public WilksLambdaTest(RealMatrix hypothesisSumOfSquares,
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
        this.lambda = getWilksLambda();

        this.totalSampleSize = totalSampleSize;
        // calculate nuE
        this.nuE = totalSampleSize - designRank;
    }
    
    
    /**
     * Calculate the denominator degrees of freedom for the WL, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getDenominatorDF()
    {                
        double df = Double.NaN;
        
        if (a*a*b*b <= 4)
        {
            df = nuE - b + 1;
        }
        else
        {
            double gDenominator = (a*a + b*b - 5);
            if (gDenominator == 0)
                throw new IllegalArgumentException("Within and between subject contrasts " +
                		"yielded divide by zero: row of C=" + a + ", cols of U=" + b);
            double g = Math.sqrt((a*a*b*b - 4) / gDenominator);
            df = (g*(nuE - (b - a +1)/2)) - (a*b - 2)/2;
        }
        
        return df;
    }

    /**
     * Calculate the non-centrality parameter for the WL, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    private double getNonCentrality()
    {
        double adjustedW = Double.NaN;
        double g = Double.NaN;

        if (a*a*b*b <= 4) 
        {
            g = 1;
            adjustedW = lambda;
        }
        else
        {
            g = Math.sqrt((a*a*b*b - 4) / (a*a + b*b - 5));
            adjustedW = Math.pow(lambda, 1/g);
        }
        
        double omega;
        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.RAO_TWO_MOMENT_OMEGA_MULT)
        {
            omega = totalSampleSize * g * (1 - adjustedW) / adjustedW;
        }
        else
        {
            omega = getDenominatorDF() * (1 - adjustedW) / adjustedW;
        }
        if (Math.abs(omega) < TOLERANCE) omega = 0;
        return Math.abs(omega);
    }

    /**
     * Calculate the numerator degrees of freedom for the WL, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return numerator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getNumeratorDF()
    {
        return a*b;
    }

    /**
     * Calculate the observed F for the WL, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return observed F
     * @throws IllegalArgumentException
     */
    private double getObservedF()
    {       
        // TODO: LAMBA is different based on data or power analysis
        double association = 0.0;
        
        if (a*a*b*b <= 4) 
        {
            association = 1 - lambda;
        }
        else
        {
            double g = Math.sqrt((a*a*b*b - 4) / (a*a + b*b - 5));
            association = 1 - Math.pow(lambda, 1/g);
        }
        
        double ddf = getDenominatorDF();
        double ndf = getNumeratorDF();
        return ((association) / ndf) / ((1 - association) / ddf);
    }
    
    /**
     * Compute a Wilks Lamba statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getWilksLambda(DistributionType type)
    throws IllegalArgumentException
    {
        
        RealMatrix adjustedH = hypothesisSumOfSquares;
        if (type != DistributionType.DATA_ANALYSIS_NULL && ((s == 1 && p > 1) ||
                fMethod == FApproximation.RAO_TWO_MOMENT_OMEGA_MULT))
        {
            adjustedH = hypothesisSumOfSquares.scalarMultiply(nuE/totalSampleSize);
        }
        
        RealMatrix T = adjustedH.add(errorSumOfSquares);
        RealMatrix inverseT = new LUDecomposition(T).getSolver().getInverse();

        RealMatrix EinverseT = errorSumOfSquares.multiply(inverseT);
        
        double det = new LUDecomposition(EinverseT).getDeterminant();
        return det;
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
