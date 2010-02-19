package edu.cudenver.bios.matrix.test;

import junit.framework.TestCase;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.QRDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import edu.cudenver.bios.matrix.GramSchmidtOrthonormalization;

public class TestMatrixOrthonormalization
{
    public static void main(String args[])
    {
        double[][] data = {{1,-1},{-1,0},{0,-1}};

       // RealVector v = new ArrayRealVector(data);
        RealMatrix m = new Array2DRowRealMatrix(data);
        GramSchmidtOrthonormalization norm = new GramSchmidtOrthonormalization(m);
        RealMatrix Q = norm.getQ();
        System.out.println(Q);
        System.out.println(norm.getR());
//
        System.out.println(Q.transpose().multiply(Q));
        System.out.println(Q.multiply(norm.getR()));
        
        System.out.println(new QRDecompositionImpl(m).getQ());

        double norm2 = m.getFrobeniusNorm();
        System.out.println("Norm=" + norm2);
    }
}
