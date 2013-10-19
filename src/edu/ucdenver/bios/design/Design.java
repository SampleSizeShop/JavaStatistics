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
package edu.ucdenver.bios.design;

import edu.cudenver.bios.power.parameters.PowerParameters;

/**
 * Class describing the study design for a linear model.  
 * <p>
 * The following matrices should be specified
 * <ul>
 * <li>Design matrix - (required) essence matrix for the model design</li>
 * <li>Beta matrix - (required) estimated regression coefficients matrix.  When a baseline
 * covariate is specified, this matrix will have both fixed and random components</li>
 * <p></p>
 * For designs with fixed effects only, a single covariance matrix should be specified
 * <ul>
 * <li>Sigma (error)  matrix - estimated covariance matrix for the residual error</li>
 * </ul>
 * For designs with a baseline covariate, three covariance matrices should be specified:
 * <ul>
 * <li>Sigma Outcome (Y) - covariance matrix of responses</li>
 * <li>Sigma Gaussian Random (G) - covariance matrix of the baseline covariate (1x1)</li>
 * <li>SIgma Outcome Gaussian Random (YG) - covariance between responses and baseline covariate</li>
 * </ul>
 * The user may perform multiple power calculations on the above set of
 * matrices by specifying lists (with the corresponding "add" function) of the following values:
 * <ul>
 * <li>alpha - type I error values</li>
 * <li>group sample size - size of each group</li>
 * <li>power - desired power (for sample size calculations)</li>
 * <li>test - statistical test (ex. Wilk's Lambda, Hotelling-Lawley Trace, etc.)</li>
 * <li>betaScale - scale factors for the beta matrix.  Allows tweaking of mean difference estimates</li>
 * <li>sigmaScale - scale factors for the sigma error matrix.  Allows tweaking of residual covariance estimates</li>
 * <li>Power method - for designs with a baseline covariate, the user request unconditional or quantile power</li>
 * <li>Quantile - for quantile power calculations, the user may enter various quantiles (ex. 0.5 = median power)</li>
 * </ul>
 * @see PowerParameters
 * @see edu.cudenver.bios.matrix.DesignEssenceMatrix
 * @see edu.cudenver.bios.matrix.FixedRandomMatrix
 * @author Sarah Kreidler
 *
 */
public class Design {
    
}
