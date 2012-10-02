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
package edu.cudenver.bios.matrix.test;

import junit.framework.TestCase;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Precision;

import edu.cudenver.bios.matrix.GramSchmidtOrthonormalization;

/**
 * Unit test for GramSchmidtOrthonormalization.  Verifies the following
 * identities for the produced Q,R matrices:
 * <ul>
 * <li>A = QxR</li>
 * <li>Q'Q = I</li>
 *  </ul>
 * 
 * @see edu.cudenver.bios.matrix.GramSchmidtOrthonormalization
 * @author Sarah Kreidler
 *
 */
public class TestMatrixOrthonormalization extends TestCase
{
    protected static final double TOLERANCE = 0.00001; 
    protected static final double[][] data = {{1,2},{-4,0},{0,2}};
    // create a generic matrix
    protected RealMatrix A = new Array2DRowRealMatrix(data);
    // orthonormalize it
    protected GramSchmidtOrthonormalization norm = new GramSchmidtOrthonormalization(A);
    
    /**
     * Verify that the A = QxR for the Q, R matrices produced by the
     * orthonormalization
     */
    public void testQRshouldBeA()
    {
        RealMatrix Q = norm.getQ();
        RealMatrix R = norm.getR();

        // verify that QR = A
        RealMatrix shouldBeOriginalMatrix = Q.multiply(R);
        // Make sure that QR has same dimension as original A matrix
        if (shouldBeOriginalMatrix.getRowDimension() != data.length ||
                shouldBeOriginalMatrix.getColumnDimension() != data[0].length)
        {
            fail();
        }
        // make sure elements of QR match elements of A (within tolerance)
        for(int r = 0; r < shouldBeOriginalMatrix.getRowDimension(); r++)
        {
            for(int c = 0; c < shouldBeOriginalMatrix.getColumnDimension(); c++)
            {
                if (Precision.compareTo(shouldBeOriginalMatrix.getEntry(r, c), data[r][c], TOLERANCE) != 0) 
                    fail();
            }
        }
        
        assertTrue(true);
    }
    
    /**
     * Verify that the Q'Q = I for the Q matrix produced by the
     * orthonormalization
     */
    public void testQQisIdentity()
    {
        RealMatrix Q = norm.getQ();
        
        // verify that Q'Q = identity
        RealMatrix shouldBeIdentityMatrix = Q.transpose().multiply(Q);
        // make sure the matrix is sqaure
        if (!shouldBeIdentityMatrix.isSquare())
        {
            fail();
        }
        // make sure the diagonal elements are one (within tolerance), and off diagonals
        // are zero (within tolerance)
        for(int r = 0; r < shouldBeIdentityMatrix.getRowDimension(); r++)
        {
            for(int c = 0; c < shouldBeIdentityMatrix.getColumnDimension(); c++)
            {
                double shouldBeValue = (r == c) ? 1 : 0;
                if (Precision.compareTo(shouldBeIdentityMatrix.getEntry(r, c), shouldBeValue, TOLERANCE) != 0) 
                    fail();
            }
        }
        assertTrue(true);
    }
}
