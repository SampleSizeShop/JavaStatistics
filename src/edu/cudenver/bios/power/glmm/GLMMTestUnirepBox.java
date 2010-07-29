package edu.cudenver.bios.power.glmm;

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMTestUnirepBox extends GLMMTestUnivariateRepeatedMeasures
{
    double unirepEpsilonExpectedValue = Double.NaN;
    
    public GLMMTestUnirepBox(GLMMPowerParameters params)
    {
        // unirep base class will calculate epsilon for box correction
        super(params);
    }
    
    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix U = params.getWithinSubjectContrast();
        
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();
        
        double df = Double.NaN;

        // for the conservative, or "Box" test, we adjust by epsilon only for
        // power analysis under the alternative.  The ddf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = b*(N - r) * this.unirepEpsilon;
        else
            df = (N - r);
        
        return df;
    }

    @Override
    public double getNonCentrality(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*unirepEpsilon;
    }

    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();

        double df = Double.NaN;

        // for the conservative, or "Box" test, we adjust by epsilon only for
        // power analysis under the alternative.  The ndf are the same for power
        // under the null and for data analysis
        if (type == DistributionType.POWER_ALTERNATIVE)
            df = a * b * this.unirepEpsilon;
        else
            df = a;

        return df;
    }

}
