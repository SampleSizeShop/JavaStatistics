/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2010 Regents of the University of Colorado.  
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package edu.ucdenver.bios.statisticaltest;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class HotellingLawleyTraceTest extends StatisticalTest {

    // hypothesis sum of squares
    private RealMatrix hypothesisSumOfSquares;
    // error sum of squares
    private RealMatrix errorSumOfSquares;
    // value of the HLT statistic
    private double HLT;

    // number of rows in the between participant contrast
    private double a;
    // number of columns in the within participant contrast
    private double b;
    // total number of outcomes (i.e. columns of beta)
    private int p;
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
    public HotellingLawleyTraceTest(RealMatrix hypothesisSumOfSquares,
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

        // calculate the value of the statistic
        this.HLT = getHotellingLawleyTrace();

        // calculate nuE
        this.nuE = totalSampleSize - designRank;
    }


    /**
     * Calculate the denominator degrees of freedom for the HLT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return denominator degrees of freedom
     * @throws IllegalArgumentException
     */
    private double getDenominatorDF()
    {       
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  

        double df = Double.NaN;
        if (fMethod == FApproximation.PILLAI_ONE_MOMENT ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT)
        {
            df = s * (nuE - b -1) + 2;
        }
        else
        {
            double t1 = nuE * nuE - nuE * (2 * b + 3) + b * (b + 3);
            double t2 = (nuE * (a  + b + 1) - (a + 2 * b + b * b - 1));
            df = 4 + (a * b + 2) * (t1/t2);
        }
        // TODO Auto-generated method stub
        return df;
    }

    /**
     * Calculate the non-centrality parameter for the HLT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return non-centrality parameter
     * @throws IllegalArgumentException
     */
    private double getNonCentrality()
    {               
        // minimum of a and b dimensions
        double s = (a < b) ? a : b;  

        if ((s == 1 && p > 1) ||
                fMethod == FApproximation.PILLAI_ONE_MOMENT_OMEGA_MULT ||
                fMethod == FApproximation.MCKEON_TWO_MOMENT_OMEGA_MULT)
        {
            return nuE * HLT;
        }
        else
        {
            return getDenominatorDF() * HLT / s;
        }
    }

    /**
     * Calculate the numerator degrees of freedom for the HLT.  
     * 
     * @return numerator degrees of freedom
     */
    private double getNumeratorDF()
    {
        return a * b;
    }

    /**
     * Calculate the observed F for the HLT, based on
     * whether the null or alternative hypothesis is assumed true.  
     * 
     * @param type distribution type
     * @return observed F
     */
    public double getDataAnalysisObservedF()
    {       
        double ddf = getDenominatorDF();
        double ndf = getNumeratorDF();
        return HLT * ((nuE-b-1)*ddf) / (ndf*(ddf-2));
    }

    public double getPowerObservedF() {
        return getNonCentrality() / getNumeratorDF();
    }

    /**
     * Compute a Hotelling-Lawley Trace statistic
     * 
     * @param H hypothesis sum of squares matrix
     * @param E error sum of squares matrix
     * @returns F statistic
     */
    private double getHotellingLawleyTrace()
            throws IllegalArgumentException {
        RealMatrix inverseE = new LUDecomposition(this.errorSumOfSquares).getSolver().getInverse();
        RealMatrix HinverseE = this.hypothesisSumOfSquares.multiply(inverseE);

        return HinverseE.getTrace();
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
