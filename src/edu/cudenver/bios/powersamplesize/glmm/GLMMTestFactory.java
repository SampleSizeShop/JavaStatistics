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
            return new GLMMTestWilksLambda(params);
        case HOTELLING_LAWLEY_TRACE:
            return new GLMMTestHotellingLawley(params);
        case UNIREP:
            return new GLMMTestUnivariateRepeatedMeasures(params);
        case PILLAI_BARTLETT_TRACE:
            return new GLMMTestPillaiBartlett(params);
        default:
            throw new IllegalArgumentException("Unknown GLMM test statistic");
        }
    }
}
