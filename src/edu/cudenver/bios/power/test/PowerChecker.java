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
package edu.cudenver.bios.power.test;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.w3c.dom.Document;
import org.xml.sax.InputSource;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.utils.ConfidenceInterval;

/**
 * Helper class which runs power calculations for the java unit tests
 * and compares against simulated values and SAS output
 * @author Sarah Kreidler
 */
public class PowerChecker
{
    private class SASPower extends GLMMPower
    {
        double delta = Double.NaN;
        public SASPower(Test test, double alpha,
                double nominalPower, double actualPower, int sampleSize,
                double betaScale, double sigmaScale,
                GLMMPowerParameters.PowerMethod method,
                double quantile, ConfidenceInterval confidenceInterval, double delta)
        {
            super(test, alpha, nominalPower, actualPower, sampleSize,
                    betaScale, sigmaScale, method, quantile, confidenceInterval);
            this.delta = delta;
        }
    }

    private boolean verbose = false;
    private double positivityThreshold = CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
    private double symmetryThreshold = CholeskyDecomposition.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;
    private static final int SIMULATION_SIZE = 10000;

    private boolean simulate = true;
    private List<GLMMPower> sasPowers = null;

    // maximum allowed deviation
    private double sasTolerance = 0.00001;
    private double simTolerance = 0.05;

    // results
    private ArrayList<Result> checkerResults = new ArrayList<Result>();
    private double maxSasDeviation = 0;
    private double maxSimDeviation = 0;
    private double maxSaslowerCIDeviation = -1;
    private double maxSasUpperCIDeviation = -1;
    private int sasResultsIndex = 0;

    /**
     * class to handle timing results
     */
    private Timer timer = new Timer();
    public class Timer
    {
        public long calculationMilliseconds;
        public long simulationMilliseconds;
        public long simulationMeanMilliseconds;

        public Timer()
        {
            reset();
        };

        public void reset()
        {
            calculationMilliseconds = 0;
            simulationMilliseconds = 0;
            simulationMeanMilliseconds = 0;
        }

        public void addCalculationTime(long time)
        {
            calculationMilliseconds += time;
        }

        public void addSimulationTime(long time)
        {
            simulationMilliseconds += time;
        }

        public void addMeanSimulationTime(long time) { simulationMeanMilliseconds += time; }
    }

    /**
     * Class describing the results of power comparisons
     * againsts SAS and simulation
     */
    public class Result
    {
        GLMMPower calculatedPower;
        double sasPower;
        double sasCILower;
        double sasCIUpper;
        double simulatedPower;
        public double sasDeviation;
        double sasCILowerDeviation;
        double sasCIUpperDeviation;
        public double simulationDeviation;

        public Result(GLMMPower calculatedPower,
                double sasPower,
                double sasCILower,
                double sasCIUpper,
                double sasDeviation,
                double sasCILowerDeviation,
                double sasCIUpperDeviation,
                double simulatedPower,
                double simulationDeviation)
        {
            this.calculatedPower = calculatedPower;
            this.sasPower = sasPower;
            this.sasCILower = sasCILower;
            this.sasCIUpper = sasCIUpper;
            this.sasDeviation = sasDeviation;
            this.sasCILowerDeviation = sasCILowerDeviation;
            this.sasCIUpperDeviation = sasCIUpperDeviation;
            this.simulatedPower = simulatedPower;
            this.simulationDeviation = simulationDeviation;
        }

        @Override
        public String toString() {
            return "Result{" +
                    "calculatedPower=" + calculatedPower.toString() +
                    ", sasPower=" + sasPower +
                    ", sasCILower=" + sasCILower +
                    ", sasCIUpper=" + sasCIUpper +
                    ", simulatedPower=" + simulatedPower +
                    ", sasDeviation=" + sasDeviation +
                    ", sasCILowerDeviation=" + sasCILowerDeviation +
                    ", sasCIUpperDeviation=" + sasCIUpperDeviation +
                    ", simulationDeviation=" + simulationDeviation +
                    '}';
        }
    };

    /**
     * Constructor.
     */
    public PowerChecker() {}

    /**
     * Create a power checker which optionally does
     * not perform simulation.
     * @param compareAgainstSimulation
     */
    public PowerChecker(boolean compareAgainstSimulation)
    {
        this.simulate = compareAgainstSimulation;
    }

    /**
     * Create a power checker which compares the power
     * values against an XML formatted input file of SAS results.
     * The user may optionally disable simulation comparisons.
     *
     * @param sasPowers SAS power values
     * @param simulate if true, power values
     * are compared against simulation.
     */
    public PowerChecker(List<GLMMPower> sasPowers, boolean simulate) {
        this.sasPowers = sasPowers;
        this.simulate = simulate;
    }

    /**
     * Create a power checker which compares the power
     * values against an XML formatted input file of SAS results.
     * The user may optionally disable simulation comparisons.
     *
     *  @param sasOutputFile XML file of SAS power values
     * @param compareAgainstSimulation if true, power values
     * are compared against simulation.
     */
    public PowerChecker(String sasOutputFile, boolean compareAgainstSimulation)
    throws IllegalArgumentException
    {
        this.simulate = compareAgainstSimulation;
        this.simulate = true;
        FileReader reader = null;
        // parse the sas xml file
        try
        {
            DocumentBuilderFactory factory =  DocumentBuilderFactory.newInstance();

            DocumentBuilder builder = factory.newDocumentBuilder();
            reader = new FileReader(sasOutputFile);
            Document doc = builder.parse(new InputSource(reader));
            sasPowers = SASOutputParser.parsePowerResults(doc);
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException("Parsing of SAS XML failed: " + e.getMessage());
        }
        finally
        {
            if (reader != null) try {reader.close();} catch (Exception e) {};
        }
    }

    /**
     * Run the power calculations for the specified parameters and tests
     * and assert whether they match simulation
     *
     * @param params
     * @return the number of powers that failed to match simulation
     */
    public void checkPower(GLMMPowerParameters params)
    {
        // create a power calculator
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        calc.setPositivityThreshold(positivityThreshold);
        calc.setSymmetryThreshold(symmetryThreshold);
        // perform the calculations
        if (verbose) System.out.println("Calculating power...");
        long startTime = System.currentTimeMillis();
        List<Power> results = null;
        try
        {
            results = calc.getPower(params);
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        long calcTime = System.currentTimeMillis() - startTime;
        if (verbose) System.out.println("Done.  Elapsed time: " +  ((double) calcTime / (double) 1000) + " seconds");
        timer.addCalculationTime(calcTime);
        int reps = 10000;
        // perform the simulation if requested
        List<Power> simResults = null;
        if (simulate)
        {
            long[] simList = new long[reps];
            for (int i = 0; i < reps; i++) {
                try {
                    if (verbose) System.out.println("Simulating power...");
                    startTime = System.currentTimeMillis();

                    simResults = calc.getSimulatedPower(params, SIMULATION_SIZE);
                    long simTime = System.currentTimeMillis() - startTime;
                    if (verbose) System.out.println("Done.  Elapsed time: " +  ((double) simTime / (double) 1000) + " seconds");
                    simList[i] = simTime;

                } catch (PowerException e) {
                    System.out.println("Simulation failed: " + e.getMessage());
                }
            }
            long sum = 0;
            for (int i = 0; i < reps; i++)
                sum = sum + simList[i];

            timer.addMeanSimulationTime(sum/(long)reps);
            timer.addSimulationTime(sum);

        }

        // simResults.forEach((result) -> System.out.println(result.toString()));


        // accumulate results and calculate maximum absolute deviation
        for(int i = 0; i < results.size(); i++, sasResultsIndex++)
        {
            GLMMPower power = (GLMMPower) results.get(i);
            GLMMPower simPower = (simResults != null ? (GLMMPower) simResults.get(i) : null);
            GLMMPower sasPower = (sasPowers != null ? (GLMMPower) sasPowers.get(sasResultsIndex) : null);
            double simDeviation = Double.NaN;
            double sasDeviation = Double.NaN;
            double sasLowerCIDeviation = Double.NaN;
            double sasUpperCIDeviation = Double.NaN;

            ConfidenceInterval sasCI = null;
            ConfidenceInterval calcCI = null;

            if (simPower != null)
            {
                simDeviation = Math.abs(power.getActualPower() - simPower.getActualPower());
                if (simDeviation > maxSimDeviation) maxSimDeviation = simDeviation;
            }
            if (sasPower != null)
            {
                sasDeviation = Math.abs(power.getActualPower() - sasPower.getActualPower());
                if (sasDeviation > maxSasDeviation) maxSasDeviation = sasDeviation;
                sasCI = sasPower.getConfidenceInterval();
                calcCI = power.getConfidenceInterval();
                if (sasCI != null && calcCI != null)
                {
                    sasLowerCIDeviation = Math.abs(calcCI.getLowerLimit() - sasCI.getLowerLimit());
                    sasUpperCIDeviation = Math.abs(calcCI.getUpperLimit() - sasCI.getUpperLimit());
                    if (sasLowerCIDeviation > maxSaslowerCIDeviation) maxSaslowerCIDeviation = sasLowerCIDeviation;
                    if (sasUpperCIDeviation > maxSasUpperCIDeviation) maxSasUpperCIDeviation = sasUpperCIDeviation;
                }
            }

            checkerResults.add(new Result(power,
                (sasPower != null ? sasPower.getActualPower(): Double.NaN),
                (sasCI != null ? sasCI.getLowerLimit(): Double.NaN),
                (sasCI != null ? sasCI.getUpperLimit(): Double.NaN),
                sasDeviation,
                sasLowerCIDeviation,
                sasUpperCIDeviation,
                (simPower != null ? simPower.getActualPower() : Double.NaN),
                simDeviation));
        }
    }

    /**
     * Get the timing results
     * @return timer object
     */
    public Timer getTiming() {
        return timer;
    }

    /**
     * Get the power comparison results
     * @return list of results
     */
    public List<Result> getResults() {
        return checkerResults;
    }

    /**
     * Enable/disable verbose output
     */
    public void setVerbose(boolean enabled)
    {
        verbose = enabled;
    }

    /**
     * Clear the current comparison.  Allows reuse of the
     * same power checker object for multiple runs.
     */
    public void reset()
    {
        checkerResults.clear();
        timer.reset();
        maxSasDeviation = 0;
        maxSimDeviation = 0;
        maxSaslowerCIDeviation = -1;
        maxSasUpperCIDeviation = -1;
        sasResultsIndex = 0;
    }

    /**
     * Returns true if the max deviation from SAS
     * power values is lower than the input tolerance level.
     * This indicates that the power results match within
     * the desired tolerance.
     *
     * @param tolerance tolerance to compare to max deviation
     * @return true if power difference is within the desired
     * tolerance.
     */
    public boolean isSASDeviationBelowTolerance(double tolerance) {
        return (maxSasDeviation <= tolerance);
    }

    /**
     * Returns true if the max deviation from SAS
     * power values is lower than the tolerance level.
     * This indicates that the power results match within
     * the desired tolerance.
     * @return true if power difference is within the desired
     * tolerance.
     */
    public boolean isSASDeviationBelowTolerance()
    {
        return (maxSasDeviation <= sasTolerance);
    }

    /**
     * Returns true if the max deviation from simulated
     * power values is lower than the tolerance level.
     * This indicates that the power results match within
     * the desired tolerance.
     * @return true if power difference is within the desired
     * tolerance.
     */
    public boolean isSimulationDeviationBelowTolerance()
    {
        return (maxSimDeviation <= simTolerance);
    }

    /**
     * Set the tolerance for comparison between SAS
     * power values and calculated power values.
     * @param tolerance max allowed difference
     */
    public void setSASTolerance(double tolerance)
    {
        sasTolerance = tolerance;
    }

    /**
     * Set the tolerance for comparison between simulated
     * power values and calculated power values.
     * @param tolerance max allowed difference
     */
    public void setSimulationTolerance(double tolerance)
    {
        simTolerance = tolerance;
    }

    /**
     * Set the minimum number which is considered "positive"
     * to control for numerical instability.
     * @param threshold minimum positive number.
     */
    public void setPositivityThreshold(double threshold)
    {
        this.positivityThreshold = threshold;
    }

    /**
     * Set the tolerance for comparing matrix values which
     * should be symmetric.  Controls for numerical
     * instability in some matrix computations.
     * @param threshold tolerance for comparing symmetric
     * matrix cells.
     */
    public void setSymmetryThreshold(double threshold)
    {
        this.symmetryThreshold = threshold;
    }

    /**
     * Get the max absolute deviation from SAS
     * power values.
     * @return max absolute deviation from SAS
     */
    public double getMaxSasDeviation() {
        return maxSasDeviation;
    }

    /**
     * Get the max absolute deviation from simulated
     * power values.
     * @return max absolute deviation from simulated
     */
    public double getMaxSimDeviation() {
        return maxSimDeviation;
    }

    /**
     * Get the max absolute deviation from SAS
     * lower CI limits for power.
     * @return max absolute deviation from SAS CI lower limit.
     */
    public double getMaxSaslowerCIDeviation() {
        return maxSaslowerCIDeviation;
    }

    /**
     * Get the max absolute deviation from SAS
     * upper CI limits for power.
     * @return max absolute deviation from SAS CI upper limit.
     */
    public double getMaxSasUpperCIDeviation() {
        return maxSasUpperCIDeviation;
    }


}
