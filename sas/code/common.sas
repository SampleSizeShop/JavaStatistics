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
/* 
* Common definitions for SAS test cases.  You will need to edit this file
* to point to your local installations of POWERLIB, LINMOD, and Glueck&Muller 2003
* code in order to run these modules.
*/
/* Directory in which to write output from SAS test cases - DO NOT EDIT */
%LET OUTPUT_DATA_DIRECTORY = ..\..\data;
/* Directory containing input data for the SAS test cases - DO NOT EDIT */
%LET INPUT_DATA_DIRECTORY = ..\data;
/* Sas Test Code directory - DO NOT EDIT */
%LET CODE_DIRECTORY = ..\code;
/* Java results directory - DO NOT EDIT */
%LET RESULTS_DIRECTORY = ..\..\text\results;

/* 
* The software modules listed below were last available from 
* http://www.ehpr.ufl.edu/muller/software_agreement
* The required software includes POWERLIB and LINMOD 
*/
%LET POWERLIB_IML_FILE = C:\KeithMullerSoftware\power\Iml\POWERLIB21.IML;
%LET LINMOD_IML_DIRECTORY = C:\KeithMullerSoftware\linmod33\SOURCE;
/* The code listed below is not publicly available.  Please contact TBD to request a copy */
%LET GLUECK_MULLER_IML_DIRECTORY = C:\KeithMullerSoftware\BASE01\kem\research\randomx\baseline\iml;

