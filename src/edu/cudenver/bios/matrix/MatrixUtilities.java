/*
 * Java Statistics.  A java library providing power/sample size estimation for
 * the general linear model.
 *
 * Copyright (C) 2016 Regents of the University of Colorado.
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

import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * A class containing some matrix utility methods.
 *
 * <p>
 * These methods are slowly being migrated from the class MatrixUtils,
 * to deal more simply with the fact that there is also an Apache
 * Commons Math class called MatrixUtils.
 */
public class MatrixUtilities {
    /**
     * Construct an instance of this class.
     *
     * <p>
     * This constructor is private, to prevent instantiation by clients.
     */
    private MatrixUtilities() {
    }

    /**
     * Print a RealMatrix to standard out.
     *
     * @param label A label to print first.
     * @param rm    The RealMatrix.
     */
    public static void dump(String label, RealMatrix rm) {
        if (rm == null) {
            System.out.println("== " + label + " (null)");
            return;
        }

        System.out.println("== " + label + " (" + rm.getRowDimension() + " x " + rm.getColumnDimension() + ")");

        for (int i = 0, m = rm.getRowDimension(); i < m; ++ i) {
            System.out.print("  ");

            for (int j = 0, n = rm.getColumnDimension(); j < n; ++ j) {
                System.out.print(rm.getEntry(i, j));
                System.out.print(' ');
            }

            System.out.println();
        }

        System.out.println();
    }

    /**
     * Force a square RealMatrix to be symmetric.
     *
     * @param rm The RealMatrix.
     *
     * @return The same RealMatrix, modified if necessary
     *         to be symmetric.
     *
     * @throws NonSquareMatrixException if the RealMatrix is
     *                                  not square.
     */
    public static RealMatrix forceSymmetric(RealMatrix rm) {
        int m = rm.getRowDimension();
        int n = rm.getColumnDimension();

        if (m != n) {
            throw new NonSquareMatrixException(m, n);
        }

        for (int i = 0; i < m; ++ i) {
            for (int j = i + 1; j < n; ++ j) {
                double value = (rm.getEntry(i, j) + rm.getEntry(j, i)) / 2;
                rm.setEntry(i, j, value);
                rm.setEntry(j, i, value);
            }
        }

        return rm;
    }
}
