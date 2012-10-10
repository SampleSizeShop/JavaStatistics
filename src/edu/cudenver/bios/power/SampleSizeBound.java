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
 * A pojo describing a sample size boundary for
 * searching.
 * @author Sarah Kreidler
 *
 */
public class SampleSizeBound {
    
    public enum SampleSizeError {
        MAX_SAMPLE_SIZE_EXCEEDED,
        SAMPLE_SIZE_UNDEFINED
    }
    
    // the sample size
    int sampleSize;
    // the actual power associated with the bound
    double actualPower;
    // error code
    SampleSizeError error = null;
    
    /**
     * Constructor for successful boundaries
     * @param sampleSize
     * @param actualPower
     */
    public SampleSizeBound(int sampleSize, double actualPower) {
        this.sampleSize = sampleSize;
        this.actualPower = actualPower;
    }
    
    /**
     * Constructor indicating that sample size failed
     * @param sampleSize
     * @param actualPower
     * @param error
     */
    public SampleSizeBound(int sampleSize, double actualPower,
            SampleSizeError error) {
        this.sampleSize = sampleSize;
        this.actualPower = actualPower;
        this.error = error;
    }

    /**
     * Get the sample size value
     * @return sample size
     */
    public int getSampleSize() {
        return sampleSize;
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
    public SampleSizeError getError() {
        return error;
    }   
}
