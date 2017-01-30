/*
 * Java Statistics.  A java library providing power/sample size estimation for
 * the general linear model.
 *
 * Copyright (C) 2017 Regents of the University of Colorado.
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

import edu.cudenver.bios.utils.Supplier; // in Java 7
//import java.util.function.Supplier; // in Java 8
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
     * The line separator.
     */
    private static final String EOL = System.getProperty("line.separator");

    /**
     * Construct an instance of this class.
     *
     * <p>
     * This constructor is private, to prevent instantiation by clients.
     */
    private MatrixUtilities() {
    }

    /**
     * Create a string representing the content of a RealMatrix,
     * suitable for logging.
     *
     * @param label An optional label.
     * @param rm    The RealMatrix.
     *
     * @return The string.
     */
    public static String logMessage(String label, RealMatrix rm) {
        StringBuilder sb = new StringBuilder();

        sb.append(EOL);

        if (label != null) {
            sb.append("  ")
              .append(label)
              .append(EOL);
        }

        if (rm == null) {
            sb.append("   (null)");
            sb.append(EOL);
        } else {
            int m = rm.getRowDimension();
            int n = rm.getColumnDimension();

            sb.append("   (")
              .append(m)
              .append(" x ")
              .append(n)
              .append(")")
              .append(EOL);

            for (int i = 0; i < m; ++ i) {
                sb.append("    ");

                for (int j = 0; j < n; ++ j) {
                    sb.append(rm.getEntry(i, j))
                      .append(' ');
                }

                sb.append(EOL);
            }
        }

        return sb.toString();
    }

    /**
     * Given a label and a matrix, return a supplier of a
     * log message for that label and matrix.
     *
     * @param label The label.
     * @param rm    The matrix.
     *
     * @return The supplier.
     */
    public static Supplier<Object> logMessageSupplier(final String label, final RealMatrix rm) {
        return new Supplier<Object>() {
            @Override
            public Object get() {
                return logMessage(label, rm);
            }
        };
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
