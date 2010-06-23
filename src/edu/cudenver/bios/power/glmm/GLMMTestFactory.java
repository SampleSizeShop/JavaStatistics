package edu.cudenver.bios.power.glmm;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMTestFactory
{
    public static GLMMTest createGLMMTest(GLMMPowerParameters params)
    throws IllegalArgumentException
    {
        switch (params.getTest())
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
