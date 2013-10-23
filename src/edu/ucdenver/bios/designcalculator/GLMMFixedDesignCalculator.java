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

import edu.ucdenver.bios.criteria.PrecisionCriteria;
import edu.ucdenver.bios.criteria.RejectionCriteria;
import edu.ucdenver.bios.design.Design;
import edu.ucdenver.bios.design.GLMMFixedDesign;
import edu.ucdenver.bios.designcalculator.exception.InvalidCriteriaException;
import edu.ucdenver.bios.designcalculator.exception.InvalidDesignException;
import edu.ucdenver.bios.designcalculator.exception.MismatchedDesignCriteriaException;

public class GLMMFixedDesignCalculator extends DesignCalculator {

    private void validateDesign(GLMMFixedDesign design) 
    throws InvalidDesignException {
            if (design == null) throw new InvalidDesignException();
            
    } 
    
    @Override
    public Power getPower(Design design, RejectionCriteria criteria) 
    throws InvalidCriteriaException, InvalidDesignException, 
    MismatchedDesignCriteriaException {
        
        GLMMFixedDesign fixedDesign = (GLMMFixedDesign) design;
        
        validateDesign(fixedDesign);

        
        
        return null;
    }

    @Override
    public DetectableDifference getDetectableDifference(Design design,
            RejectionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public DetectableDifference getDetectableDifference(Design design,
            PrecisionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public DetectableDifference getDetectableDifference(Design design,
            RejectionCriteria rejectionCriteria,
            PrecisionCriteria precisionCriteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public SampleSize getSampleSize(Design design, RejectionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public SampleSize getSampleSize(Design design, PrecisionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public SampleSize getSampleSize(Design design,
            RejectionCriteria rejectionCriteria,
            PrecisionCriteria precisionCriteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public ClusterSize getClusterSize(Design design, RejectionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public ClusterSize getClusterSize(Design design, PrecisionCriteria criteria) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public ClusterSize getClusterSize(Design design,
            RejectionCriteria rejectionCriteria,
            PrecisionCriteria precisionCriteria) {
        // TODO Auto-generated method stub
        return null;
    }

}
