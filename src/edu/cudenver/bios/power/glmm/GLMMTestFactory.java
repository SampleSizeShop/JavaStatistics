package edu.cudenver.bios.power.glmm;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMTestFactory
{
    public static GLMMTest createGLMMTest(GLMMPowerParameters params)
    throws IllegalArgumentException
    {
        switch (params.getCurrentTest())
        {
        case UNIREP:
            return new GLMMTestUnivariateRepeatedMeasures(params);
        case UNIREP_BOX:
            return new GLMMTestUnirepBox(params);
        case UNIREP_GEISSER_GREENHOUSE:
            return new GLMMTestUnirepGeisserGreenhouse(params);
        case UNIREP_HUYNH_FELDT:
            return new GLMMTestUnirepHuynhFeldt(params);
        case WILKS_LAMBDA:
            return new GLMMTestWilksLambda(params);
        case HOTELLING_LAWLEY_TRACE:
            return new GLMMTestHotellingLawley(params);
        case PILLAI_BARTLETT_TRACE:
            return new GLMMTestPillaiBartlett(params);
        default:
            throw new IllegalArgumentException("Unknown GLMM test statistic");
        }
    }
}
