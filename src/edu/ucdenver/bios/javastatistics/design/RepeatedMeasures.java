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

import java.util.ArrayList;
import java.util.List;
/**
 * This is a wrapper for the spacing in repeated measures.
 * @author Uttara Sakhadeo
 *
 */
class Spacing
{
	protected List<Integer> spacingList = null;
	
	public Spacing()
	{
		spacingList = new ArrayList<Integer>();
	}		
	public Spacing(List<Integer> spacingList) {
		//super();
		this.spacingList = spacingList;
	}
	public List<Integer> getSpacingList() {
		return spacingList;
	}		
	public void setSpacingList(List<Integer> spacingList) {
		this.spacingList = spacingList;
	}		
}

/**
 * This is a wrapper for the repeated measures information.
 * @author Uttara Sakhadeo
 *
 */
public class RepeatedMeasures extends Spacing
{
	private String name = null;
	private String type = null;
	private Integer count = null;
	//private Spacing spacingList = null;
	
	public RepeatedMeasures()
	{
		//spacingList = new Spacing();
	}	
	public RepeatedMeasures(String name) 
	{
		//super();
		this.name = name;		
	}		
	public RepeatedMeasures(String name, String type, Integer count) 
	{
		//super();
		this.name = name;
		this.type = type;
		this.count = count;		
	}			
	public RepeatedMeasures(String name, String type, Integer count,
			List<Integer> spacingList) 
	{
		//super();
		this.name = name;
		this.type = type;
		this.count = count;
		this.setSpacingList(spacingList);
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public Integer getCount() {
		return count;
	}
	public void setCount(Integer count) {
		this.count = count;
	}
	public List<Integer> getSpacingList() {
		return spacingList;
	}
	public void setSpacingList(List<Integer> spacingList) {
		this.spacingList = spacingList;
	}		
}