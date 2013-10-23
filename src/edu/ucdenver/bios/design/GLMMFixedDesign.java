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
package edu.ucdenver.bios.design;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.glmm.GLMMPowerConfidenceInterval.ConfidenceIntervalType;
import edu.cudenver.bios.power.glmm.GLMMTest.FApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateCdfApproximation;
import edu.cudenver.bios.power.glmm.GLMMTest.UnivariateEpsilonApproximation;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;

public class GLMMFixedDesign extends Design {
     
    /* matrix inputs to the power calculation.  */
    // the design essence matrix 
    /* For details please see Muller & Fetterman (2002) "Regression and ANOVA" */
    RealMatrix designEssence = null;
    // caching of X'X inverse and rank since these are order(n^3) operations
    RealMatrix XtXInverse = null;
    int designRank = -1;

    // beta matrix of regression coefficients
    RealMatrix beta = null;
    // residual error matrix 
    RealMatrix sigmaError = null;
    // if true, power will be adjusted slightly to account for uncertainty in Sigma
    boolean sigmaEstimated = false;
    
    // additional parameters must be specified when Sigma is estimated
    // in order to build the confidence limits
    // lower tail probability
    double alphaLowerConfidenceLimit = -1;
    // upper tail probability
    double alphaUpperConfidenceLimit = -1;
    // 
    int sampleSizeForEstimates = 0;
    int designMatrixRankForEstimates = 0;

}
