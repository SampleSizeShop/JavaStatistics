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
package edu.ucdenver.bios.criteria;

import org.apache.commons.math3.linear.RealMatrix;

import edu.ucdenver.bios.statisticaltest.StatisticalTest;

public class RejectionCriteria {
    // the between participant contrast for fixed predictors
    RealMatrix betweenFixedContrastMatrix;
    // the between participant contrast for random predictors
    RealMatrix betweenRandomContrastMatrix;
    
    // the within participant contrast
    RealMatrix withinParticipantContrastMatrix;
    
    // the null hypothesis matrix
    RealMatrix nullHypothesisMatrix;
    
    // the alpha level 
    double alpha;
    
    // the statistical test
    StatisticalTest test;
    
    /*
     *  A subset of the following will be set depending on which value
     *  the user is solving for (power, sample size, detectable difference, etc.)
     */
    
    // per group sample size
    int perGroupSampleSize;
    
    // desired power
    double desiredPower;   
    
}
