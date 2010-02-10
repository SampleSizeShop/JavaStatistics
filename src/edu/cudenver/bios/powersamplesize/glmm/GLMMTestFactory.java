package edu.cudenver.bios.powersamplesize.glmm;

import edu.cudenver.bios.powersamplesize.parameters.LinearModelPowerSampleSizeParameters;

public class GLMMTestFactory
{
    public static GLMMTest createGLMMTest(LinearModelPowerSampleSizeParameters params)
    throws IllegalArgumentException
    {
        switch (params.getTestStatistic())
        {
        case WILKS_LAMBDA:
            return new WilksLambdaGLMMTest(params);
        case HOTELLING_LAWLEY_TRACE:
            return new HotellingLawleyGLMMTest(params);
        case UNIREP:
            return new UnivariateRepeatedMeasuresGLMMTest(params);
        case PILLAI_BARTLETT_TRACE:
            return new PillaiBartlettGLMMTest(params);
        default:
            throw new IllegalArgumentException("Unknown GLMM test statistic");
        }
    }
}
