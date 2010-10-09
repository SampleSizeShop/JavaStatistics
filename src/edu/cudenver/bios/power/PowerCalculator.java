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
package edu.cudenver.bios.power;

import java.util.List;

import edu.cudenver.bios.power.parameters.PowerParameters;

/**
 * Abstract base class for power calculators.  All calculators should 
 * provide functionality to calculate power, simulate power, and
 * determine sample size and effect size.
 * 
 * @author Sarah Kreidler
 *
 */
public interface PowerCalculator
{
	/**
	 * Calculate power
	 * @param params
	 * @return list of power results
	 */
	List<Power> getPower(PowerParameters params);
	
	/**
	 * Determine sample size
	 * @param params
	 * @return list of sample size results
	 */
	List<Power> getSampleSize(PowerParameters params);
	
	/**
	 * Determine effect size
	 * @param params
	 * @return list of effect size results
	 */
	List<Power> getDetectableDifference(PowerParameters params);
	
	/**
	 * Determine power by simulation
	 * @param params
	 * @param iterations
	 * @return
	 */
	List<Power> getSimulatedPower(PowerParameters params, int iterations);
}
