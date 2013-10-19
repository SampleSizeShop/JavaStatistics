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

/**
 * Applies a set of quality criteria to a study design description to
 * produce a complete study design
 * 
 * @author kreidles
 *
 */
public class DesignCalculator {
    
    public Design getPower(DesignDescription designDescription, CriteriaRejection criteria) {
        // TODO
        return null;
    }
    
    public Design getDetectableDifference(DesignDescription designDescription, CriteriaRejection criteria) {
        // TODO
        return null;
    }
    
    public Design getDetectableDifference(DesignDescription designDescription, CriteriaPrecision criteria) {
        // TODO
        return null;
    }
    
    
    public Design getDetectableDifference(DesignDescription designDescription, CriteriaRejection rejectionCriteria, 
            CriteriaPrecision precisionCriteria) {
        // TODO
        return null;
    }
    
    
    public Design getSampleSize(DesignDescription designDescription, CriteriaRejection criteria) {
        // TODO
        return null;
        
    }
    
    public Design getSampleSize(DesignDescription designDescription, CriteriaPrecision criteria) {
        // TODO
        return null;
        
    }
    
    public Design getSampleSize(DesignDescription designDescription, CriteriaRejection rejectionCriteria, 
            CriteriaPrecision precisionCriteria) {
        // TODO
        return null;
        
    }
    
    public Design getClusterSize(DesignDescription designDescription, CriteriaRejection criteria) {
        // TODO
        return null;
        
    }
    
    public Design getClusterSize(DesignDescription designDescription, CriteriaPrecision criteria) {
        // TODO
        return null;
        
    }
    
    public Design getClusterSize(DesignDescription designDescription, CriteriaRejection rejectionCriteria, 
            CriteriaPrecision precisionCriteria) {
        // TODO
        return null;
        
    }
}
