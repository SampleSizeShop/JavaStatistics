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
 * Container class for meta data associated with a matrix row.
 * Primarily used with a design "essence" matrix to specify the relative
 * size of the group represented by a row in the full design matrix.
 * 
 * @see DesignEssenceMatrix
 * @author Sarah Kreidler
 *
 */
public class RowMetaData
{    
    int ratio = 1;
    /**
     * Constructor, creates an empty row meta data object
     */
    public RowMetaData() {}
    
    /**
     * Constructor.  Create a RowMetaData object with the specified
     * relative group size (i.e. group size ratio)
     * 
     * @param ratio relative group size.
     * @throws IllegalArgumentException
     */
    public RowMetaData (int ratio) 
    throws IllegalArgumentException
    {
        if (ratio <= 0) 
            throw new IllegalArgumentException("Relative group size must be greater than zero");
        this.ratio = ratio;
    }

    /**
     * Get the relative group size for the current row
     * 
     * @return relative group size
     */
	public int getRatio()
	{
		return ratio;
	}

	/**
	 * Set the relative group size for the current row
	 * 
     * @param ratio relative group size.
	 */
	public void setRatio(int ratio)
	{
		this.ratio = ratio;
	}   
    
    
}
