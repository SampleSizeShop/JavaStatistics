package edu.cudenver.bios.power.glmm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;

import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMTestUnirepHuynhFeldt extends GLMMTestUnivariateRepeatedMeasures
{
    protected static final double TOLERANCE = 0.000000000001;
    double unirepEpsilon = Double.NaN;
    double unirepEpsilonExpectedValue = Double.NaN;

    public GLMMTestUnirepHuynhFeldt(GLMMPowerParameters params)
    {
        super(params);

        // verify that U is orthonormal to an identity matrix
        // if not, build an orthonormal U from the specified U matrix
        createOrthonormalU();

        // pre-calculate the values for epsilon (correction for violation of sphericity)
        calculateUnirepCorrection();
    }

    @Override
    public double getDenominatorDF(DistributionType type)
    {
        RealMatrix U = params.getWithinSubjectContrast();
        
        // b = #columns in within subject contrast matrix
        int b = U.getColumnDimension();
        
        double df = Double.NaN;

        // HF correction, we multiply the ddf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = b*(N - r)*this.unirepEpsilonExpectedValue;
        else
            df = b*(N - r)*this.unirepEpsilon;

        return df;
    }

    @Override
    public double getNonCentrality(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        // calculate non-centrality and adjust for sphericity 
        return a*b*getObservedF(type)*unirepEpsilon;
    }

    @Override
    public double getNumeratorDF(DistributionType type)
    {
        double a = params.getBetweenSubjectContrast().getCombinedMatrix().getRowDimension();
        double b = params.getWithinSubjectContrast().getColumnDimension();
        
        double df = Double.NaN;

        // HF correction, we multiply the ndf by the epsilon estimate for
        // power analysis (under alternative) and for data analysis.  For power under
        // the null, we multiply by the expected value of the epsilon estimate
        if (type == DistributionType.POWER_NULL)
            df = a * b * this.unirepEpsilonExpectedValue;
        else
            df = a * b * this.unirepEpsilon;
        
        return df;
    }

    private void calculateUnirepCorrection()
    {          
        RealMatrix U = params.getWithinSubjectContrast();
        int b = new SingularValueDecompositionImpl(U).getRank();
        // get the sigmaStar matrix: U' *sigmaError * U
        RealMatrix sigmaStar = U.transpose().multiply(params.getScaledSigmaError().multiply(U));
        // ensure symmetry
        sigmaStar = sigmaStar.add(sigmaStar.transpose()).scalarMultiply(0.5); 
        // normalize
        sigmaStar = sigmaStar.scalarMultiply(1/sigmaStar.getTrace());

        // get the eigen values of the normalized sigmaStar matrix
        double[] eigenValues = new EigenDecompositionImpl(sigmaStar, TOLERANCE).getRealEigenvalues();
        if (eigenValues.length <= 0) throw new IllegalArgumentException("Failed to compute eigenvalues for sigma* matrix");
        Arrays.sort(eigenValues); // put eigenvalues in ascending order

        // calculate epsilon (correction for violation of sphericity)
        // to avoid looping over the eigenvalues twice, we also calculate the multiplicity for distinct eigenvalues

        // list of distinct eigenvalues with multiplicity
        ArrayList<EigenValueMultiplicityPair> distinctEigenValues = new ArrayList<EigenValueMultiplicityPair>();
        // initialize values for the first eigen value
        double first = eigenValues[0];
        distinctEigenValues.add(new EigenValueMultiplicityPair(first, 1));
        double sumLambda = first;
        double sumLambdaSquared = first * first;

        // loop over remaining eigen values, saving distinct eigen values
        for(int i = 1; i < eigenValues.length; i++)
        {
            double value = eigenValues[i];
            // build the sum & sum of squares of eigen values
            sumLambda += value;
            sumLambdaSquared += value * value;

            // determine if this is a distinct eigen value and calculate multiplicity
            EigenValueMultiplicityPair prev = distinctEigenValues.get(distinctEigenValues.size()-1);
            if (Math.abs(prev.eigenValue - value) > TOLERANCE)
            {
                // found new distinct eigen value
                distinctEigenValues.add(new EigenValueMultiplicityPair(value, 1));
            }
            else
            {
                // repeat of same eigenvalue, so  increment the multiplicity
                prev.multiplicity++;
            }
        }

        // calculate estimate of epsilon (correction for violation of spehericity assumption)
        unirepEpsilon = (sumLambda*sumLambda) / (b * (sumLambdaSquared));

        // For Huynh Feldt, we also need the expected value of the epsilon
        // estimate.  Note, the estimates of epsilon represent different functions of the
        // eigenvalues, so the resulting derivatives and expected values are specific to each test

        // calculate the expected value of the epsilon estimate
        // E[h(lambda)] = h(lambda) + g1 / (N - r)
        // h(lambda) = h1(lambda) / (b*h2(lambda)
        // see Muller, Barton (1989) for details
        double h1 = N * sumLambda * sumLambda - 2 * sumLambdaSquared;
        double h2 = (N - r) * sumLambdaSquared - (sumLambda * sumLambda);
        double g1 = 0;
        for(int i = 0; i < distinctEigenValues.size(); i++)
        {
            EigenValueMultiplicityPair evmI = distinctEigenValues.get(i);
            // derivatives of sub-equations comprising epsilon estimator
            double h1firstDerivative = (2 * N * sumLambda) - (4 *evmI.eigenValue); 
            double h1secondDerivative = 2 * N - 4;

            double h2firstDerivative = (2 * (N - r) * evmI.eigenValue) - (2 * sumLambda); 
            double h2secondDerivative = 2 * (N - r) - 2;

            // derivatives of estimate of epsilon
            double firstDerivative = ((h1firstDerivative) - ((h1 * h2firstDerivative) / h2)) / (h2 *b); 
            double secondDerivative = 
                (h1secondDerivative - ((2*h1firstDerivative*h2firstDerivative)/(h2)) +  
                        (2 * h1 * h2firstDerivative * h2firstDerivative)/(h2*h2) - 
                        (h1*h2secondDerivative)/(h2)) / 
                        (h2*b);

            // accumulate the first term of g1 (sum over distinct eigen vals of 1st derivative * eigen val ^2 * multiplicity)
            g1 += secondDerivative * evmI.eigenValue * evmI.eigenValue * evmI.multiplicity;
            // loop over elements not equal to current
            for(int j = 0; j < distinctEigenValues.size(); j++)
            {
                if (i != j)
                {
                    EigenValueMultiplicityPair evmJ = distinctEigenValues.get(j);
                    // accumulate second term of g1
                    g1 += ((firstDerivative * evmI.eigenValue * evmI.multiplicity * evmJ.eigenValue * evmJ.multiplicity) /
                            (evmI.eigenValue - evmJ.eigenValue));
                }
            }
        }

        this.unirepEpsilonExpectedValue = (N*b*unirepEpsilon - 2)/(b*(N - r - b*unirepEpsilon))  + g1 / (N - r);


        // ensure that expected value is within bounds 1/b to 1
        if (unirepEpsilonExpectedValue != Double.NaN)
        {
            if (unirepEpsilonExpectedValue < 1/b)
            {
                unirepEpsilonExpectedValue = 1/b;
            }
            else if (unirepEpsilonExpectedValue > 1)
            {
                unirepEpsilonExpectedValue = 1;
            }
        }

    }

}
