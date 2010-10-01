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
package edu.cudenver.bios.matrix;

/**
 * Container class for column specific meta data in an essence matrix.
 * Allows the user to specify  a mean and variance of the predictor 
 * distribution for a given column
 * 
 *  @see DesignEssenceMatrix
 * @author Sarah Kreidler
 *
 */
public class RandomColumnMetaData
{
	double mean = Double.NaN;
    double variance = Double.NaN;
    
    /**
     * Constructor.  Create a column meta data object with the specified predictor
     * type.  Mean and variance should be specified for random predictors.
     * 
     * @param mean mean of the predictor distribution (assumed normal)
     * @param variance variance of the predictor distribution
     */
    public RandomColumnMetaData (double mean, double variance) 
    {
        this.mean = mean;
        this.variance = variance;
    }

    /**
     * Get the mean for this column meta data type.  This value is
     * valid only for random predictors.
     * 
     * @return mean of the predictor distribution
     */
    public double getMean()
    {
        return mean;
    }

    /**
     * Set the mean for the predictor distribution.  Only valid 
     * for random predictors.
     * 
     * @param mean mean of the predictor distribution
     */
    public void setMean(double mean)
    {
        this.mean = mean;
    }

    /**
     * Get the variance of the predictor distribution.  Only valid 
     * for random predictors.
     * 
     * @return variance of predictor distribution
     */
    public double getVariance()
    {
        return variance;
    }

    /**
     * Set the variance of the predictor distribution.  Only valid for
     * random predictors.
     * 
     * @param variance variance of the predictor distribution
     */
    public void setVariance(double variance)
    {
        this.variance = variance;
    }
    
}
