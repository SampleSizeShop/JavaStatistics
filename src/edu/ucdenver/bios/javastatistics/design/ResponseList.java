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
import java.util.UUID;
/**
 * This is a wrapper for the lists of different types.
 * @author Uttara Sakhadeo
 *
 */
public class ResponseList 
{
	private UUID studyUUID;
	private String name = null;
	private List<String> dataList = null;
	
	public ResponseList()
	{
		this.dataList = new ArrayList<String>();
	}
	
	public ResponseList(String name, List<String> list)
	{
		//this.dataList = new ArrayList<String>(list);
		super();
		this.name = name;	
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public List<String> getDataList() {
		return dataList;
	}

	public void setDataList(List<String> dataList) {
		this.dataList = dataList;
	}	
	
	public String getElement(int index)
	{
		return this.dataList.get(index);
	}
	
	public void setEntry(int index, String element) {		
		this.dataList.set(index,element);
	}
	
	public String getEntry(int index) {		
		return this.dataList.get(index);
	}

	public UUID getStudyUUID() {
		return studyUUID;
	}

	public void setStudyUUID(UUID studyUUID) {
		this.studyUUID = studyUUID;
	}
	
	public void addEntry(String element) {		
		this.dataList.add(element);
	}		
}
