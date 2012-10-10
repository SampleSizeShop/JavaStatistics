/*
 * Java Statistics.  A java library providing power/sample size estimation for 
 * the general linear model.
 * 
 * Copyright (C) 2012 Regents of the University of Colorado.  
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
package edu.cudenver.bios.power;

/**
 * A pojo describing a detectable difference boundary for
 * searching.
 * @author Sarah Kreidler
 *
 */
public class DetectableDifferenceBound {
    
    public enum DetectableDifferenceError {
        MAX_BETA_SCALE_EXCEEDED,
        BETA_SCALE_UNDEFINED
    }
    
    // the beta scale representing the difference
    double betaScale = 0;
    // the actual power associated with the bound
    double actualPower = 0;
    // error code
    DetectableDifferenceError error = null;
    
    /**
     * Constructor for successful boundaries
     * @param betaScale
     * @param actualPower
     */
    public DetectableDifferenceBound(double betaScale, double actualPower) {
        this.betaScale = betaScale;
        this.actualPower = actualPower;
    }
    
    /**
     * Constructor indicating that detectable difference failed
     * @param beta scale
     * @param actualPower
     * @param error
     */
    public DetectableDifferenceBound(double betaScale, double actualPower,
            DetectableDifferenceError error) {
        this.betaScale = betaScale;
        this.actualPower = actualPower;
        this.error = error;
    }

    /**
     * Get the beta scale value
     * @return beta scale
     */
    public double getBetaScale() {
        return betaScale;
    }

    /**
     * Get the actual power
     * @return power
     */
    public double getActualPower() {
        return actualPower;
    }   
    
    /**
     * Get the actual power
     * @return power
     */
    public DetectableDifferenceError getError() {
        return error;
    }   
}
