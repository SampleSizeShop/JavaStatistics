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

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Perform column-wise Gram Schmidt orthonormalization on a
 * matrix.  The matrix A (mxn) is decomposed into two matrices 
 * Q (mxn), R (nxn) such that
 * <ul>
 * <li>A = QR
 * <li>Q'Q = Identity
 * <li>R is upper triangular
 * </ul> 
 * 
 * @author Sarah Kreidler
 *
 */
public class GramSchmidtOrthonormalization
{
    protected RealMatrix R;
    protected RealMatrix Q;

    /**
     * Perform Gram Schmidt Orthonormalization on the specified 
     * matrix. The matrix A (mxn) is decomposed into two matrices 
     * Q (mxn), R (nxn) such that
     * <ul>
     * <li>A = QR
     * <li>Q'Q = Identity
     * <li>R is upper triangular
     * </ul> 
     * The resulting Q, R matrices can be retrieved with the getQ()
     * and getR() functions.
     * 
     * @param matrix
     */
    public GramSchmidtOrthonormalization(RealMatrix matrix)
    {
        if (matrix == null) throw new IllegalArgumentException("Null matrix");
        
        // create the Q, R matrices
        int m = matrix.getRowDimension();
        int n = matrix.getColumnDimension();
        Q = MatrixUtils.createRealMatrix(m, n);
        R = MatrixUtils.createRealMatrix(n, n);

        // perform Gram Schmidt process using the following algorithm
        // let w<n> be the resulting orthonormal column vectors
        // let v<n> be the columns of the incoming matrix
        // w1 = (1/norm(v1))*v1
        // ...
        // wj = 1/norm(vj - projectionVj-1Vj)*[vj - projectionVj-1Vj]
        // where projectionVj-1Vj = (w1 * vj) * w1 + (w2 * vj) * w2 + ... + (wj-1 * vj) * wj-1
        //
        for(int i = 0; i < n; i++)
        {
            RealMatrix v = matrix.getColumnMatrix(i);
            for(int j = 0; j < i; j++)
            {
                RealMatrix Qj = Q.getColumnMatrix(j);
                double value = Qj.transpose().multiply(v).getEntry(0, 0);
                R.setEntry(j, i, value);
                v = v.subtract(Qj.scalarMultiply(value));
            }
            double norm = v.getFrobeniusNorm();
            R.setEntry(i, i, norm);
            Q.setColumnMatrix(i, v.scalarMultiply(1/norm));
        }
    }
    
    /**
     * Retrieve the orthogonal Q matrix from the orthonormalization.
     * 
     * @return Q matrix
     */
    public RealMatrix getQ()
    {
        return Q;
    }
    
    /**
     * Retrieve the upper triangular R matrix from the orthonormalization.
     * 
     * @return Q matrix
     */
    public RealMatrix getR()
    {
        return R;
    }
    
}
