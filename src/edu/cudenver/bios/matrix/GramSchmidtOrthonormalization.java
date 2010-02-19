package edu.cudenver.bios.matrix;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public class GramSchmidtOrthonormalization
{
    protected RealMatrix R;
    protected RealMatrix Q;

    public GramSchmidtOrthonormalization(RealMatrix matrix)
    {
        if (matrix == null) throw new IllegalArgumentException("Null matrix");
        
        // create the Q, R matrices
        int m = matrix.getRowDimension();
        int n = matrix.getColumnDimension();
        Q = MatrixUtils.createRealMatrix(m, n);
        R = MatrixUtils.createRealMatrix(n, n);

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
    
    public RealMatrix getQ()
    {
        return Q;
    }
    
    public RealMatrix getR()
    {
        return R;
    }
    
    
}
