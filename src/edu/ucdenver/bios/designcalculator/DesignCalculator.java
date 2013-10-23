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
package edu.ucdenver.bios.designcalculator;

import org.apache.commons.math3.linear.RealMatrix;

import edu.ucdenver.bios.criteria.PrecisionCriteria;
import edu.ucdenver.bios.criteria.RejectionCriteria;
import edu.ucdenver.bios.design.Design;
import edu.ucdenver.bios.designcalculator.exception.InvalidCriteriaException;
import edu.ucdenver.bios.designcalculator.exception.InvalidDesignException;
import edu.ucdenver.bios.designcalculator.exception.MismatchedDesignCriteriaException;

/**
 * Applies a quality criterion to a study design description to
 * produce a complete study design
 * 
 * @author kreidles
 *
 */
public abstract class DesignCalculator {
    
    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    public RealMatrix getHypothesisSumOfSquares(RealMatrix C, RealMatrix beta, RealMatrix U,
            RealMatrix thetaNull, RealMatrix XtXInverse)
    {        
        // M = C(X'X)^-1C'
        RealMatrix M = C.multiply(XtXInverse.multiply(C.transpose()));
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(beta.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
        
        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(M.multiply(thetaDiff));
        
        return hss;
        
    }
    
    /**
     * Calculate the sum of squares error matrix (the E matrix)
     * 
     * @param params matrices input by the user
     * @return error sum of squares
     */
    public RealMatrix getErrorSumOfSquares(RealMatrix U, RealMatrix sigmaError,
            double nuError)
    {        
        return U.transpose().multiply(sigmaError.multiply(U)).scalarMultiply(nuError);
    }    
    
    /**
     * 
     */
    /**
     * Calculate the power of the specified design based on the given rejection criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract Power getPower(Design design, RejectionCriteria criteria) 
    throws InvalidCriteriaException, InvalidDesignException, 
    MismatchedDesignCriteriaException;
    
    /**
     * Calculate the minimum detectable difference required to meet the specified 
     * rejection criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract DetectableDifference getDetectableDifference(Design design, RejectionCriteria criteria);
    
    /**
     * Calculate the minimum detectable difference required to meet the specified 
     * precision criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract DetectableDifference getDetectableDifference(Design design, PrecisionCriteria criteria);
    
    /**
     * Calculate the minimum detectable difference required to meet the specified 
     * rejection and precision criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract DetectableDifference getDetectableDifference(Design design, RejectionCriteria rejectionCriteria, 
            PrecisionCriteria precisionCriteria);
    
    /**
     * Calculate the minimum sample size required to meet the specified 
     * rejection criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract SampleSize getSampleSize(Design design, RejectionCriteria criteria);
    
    /**
     * Calculate the minimum sample size required to meet the specified 
     * precision criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract SampleSize getSampleSize(Design design, PrecisionCriteria criteria);
    
    /**
     * Calculate the minimum sample size required to meet the specified 
     * rejection and precision criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract SampleSize getSampleSize(Design design, RejectionCriteria rejectionCriteria, 
            PrecisionCriteria precisionCriteria);
    
    /**
     * Calculate the minimum cluster size required to meet the specified rejection criteria
     * @param design
     * @param criteria
     * @return
     */
    public abstract ClusterSize getClusterSize(Design design, RejectionCriteria criteria);
    
    /**
     * Calculate the minimum cluster size required to meet the specified precision criteria
     * 
     * @param design
     * @param criteria
     * @return
     */
    public abstract ClusterSize getClusterSize(Design design, PrecisionCriteria criteria);
    
    /**
     * Calculate the minimum cluster size required to meet the specified rejection and
     * precision criteria
     * 
     * @param design
     * @param rejectionCriteria
     * @param precisionCriteria
     * @return
     */
    public abstract ClusterSize getClusterSize(Design design, RejectionCriteria rejectionCriteria, 
            PrecisionCriteria precisionCriteria);
}
