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

import org.apache.commons.math3.linear.RealMatrix;

public class GLMMFixedRandomDesign extends Design {
    /* matrix inputs to the power calculation.  */
    // the design essence matrix 
    /* For details please see Muller & Fetterman (2002) "Regression and ANOVA" */
    private RealMatrix designEssence = null;
    // caching of X'X inverse and rank since these are order(n^3) operations
    private RealMatrix XtXInverse = null;
    private int designRank = -1;

    // beta matrix of regression coefficients associated with fixed predictors
    RealMatrix betaFixed = null;
    // beta matrix of regression coefficients associated with random predictors
    RealMatrix betaRandom = null;
    
    // residual error matrix - calculated from the remaining elements
    private RealMatrix sigmaError = null;
    // covariance for outcomes (assuming no covariates)
    private RealMatrix sigmaOutcomes = null;
    // covariance of Gaussian predictors
    private RealMatrix sigmaRandomPredictors = null;
    // covariance of outcomes and random predictors
    private RealMatrix sigmaOutcomesRandomPredictors = null;
    
    // if true, power will be adjusted slightly to account for uncertainty in Sigma
    private  boolean sigmaEstimated = false;
    
    // additional parameters must be specified when Sigma is estimated
    // in order to build the confidence limits
    // lower tail probability
    private double alphaLowerConfidenceLimit = -1;
    // upper tail probability
    private  double alphaUpperConfidenceLimit = -1;
    // 
    private int sampleSizeForEstimates = 0;
    private int designMatrixRankForEstimates = 0;
    
    
}
