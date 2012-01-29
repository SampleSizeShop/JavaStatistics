/*
 * Study Design Service for the GLIMMPSE Software System.  
 * This service stores study design definitions for users of the GLIMMSE interface.
 * Service contain all information related to a power or sample size calculation.  
 * The Study Design Service simplifies communication between different screens in the user interface.
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
package edu.ucdenver.bios.javastatistics.design;
/**
 * This is a wrapper for the clustering information.
 * @author Uttara Sakhadeo
 *
 */
public class ClusterNode 
{		
	/*--------------------
	 * Member Variables
	 *--------------------*/
	private String clusterName = null;	
	private Integer sampleSize = null;	
	/*--------------------
	 * Constructors
	 *--------------------*/	
	public ClusterNode() {}	
	/**
	 * @param clusterName
	 */
	public ClusterNode(String clusterName) {
		this.clusterName = clusterName;
	}
	/**
	 * @param clusterName
	 * @param sampleSize
	 */
	public ClusterNode(String clusterName, Integer sampleSize) {
		this.clusterName = clusterName;
		this.sampleSize = sampleSize;
	}	
	/*--------------------
	 * Getter/Setter Methods
	 *--------------------*/	
	public String getClusterName() {
		return clusterName;
	}
	public void setClusterName(String clusterName) {
		this.clusterName = clusterName;
	}	
	public Integer getSampleSize() {
		return sampleSize;
	}
	public void setSampleSize(int sampleSize) {
		this.sampleSize = new Integer(sampleSize);
	}		
}