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

/**
 * Generic Exception class for power calculations
 * @author Sarah Kreidler
 *
 */
public class PowerException extends Exception
{
	static final long serialVersionUID = -1L;
	protected PowerErrorEnum errorCode = null;
	
	/**
	 * Constructor taking an error message
	 * @param msg
	 */
	public PowerException(String msg)
	{
		super(msg);
	}
	
	/**
	 * Constructor taking an error message and error code
	 * @param msg
	 */
	public PowerException(String msg, PowerErrorEnum errorCode)
	{
	    super(msg);
	    this.errorCode = errorCode;
	}
	
	/**
	 * Constructor taking a Throwable object
	 * @param e
	 */
	public PowerException(Throwable e)
	{
		super(e);
	}

	/**
	 * Get the error code for this exception
	 * @return error code
	 */
    public PowerErrorEnum getErrorCode() {
        return errorCode;
    }	
	
}
