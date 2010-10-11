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
package edu.cudenver.bios.power.glmm;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

/**
 * Factory class for generating a GLMM test object
 * 
 * @see GLMMTest
 * @author Sarah Kreidler
 *
 */
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
