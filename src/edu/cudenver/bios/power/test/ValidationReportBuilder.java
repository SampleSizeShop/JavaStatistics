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

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;

import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker.Result;
import edu.cudenver.bios.utils.ConfidenceInterval;

/**
 * Utility routines for creating validation reports in Latex
 * and writing summary results to stdout.
 * @author Sarah Kreidler
 *
 */
public class ValidationReportBuilder {
    // version of the library.  Updated during the build
    private static final String LIBRARY_VERSION = "UNKNOWN";
    // section headers and static text
    private static final String TITLE_PREFIX = "GLIMMPSE Validation Report: ";
    private static final String SECTION_INTRO = "Introduction";
    private static final String SECTION_DESIGN = "Study Design";
    private static final String SECTION_RESULTS = "Validation Results";   
    private static final String SUBSECTION_INPUTS = "Inputs to the Power Calculation";
    private static final String SUBSECTION_TIMING = "Timing";
    private static final String SUBSECTION_SUMMARY = "Summary Statistics";
    private static final String SUBSECTION_RESULTS = "Full Validation Results";
    private static final String TEXT_INTRO1 = "The following report contains validation results " +
            "for the JavaStatistics library, a component of the GLIMMPSE software system.  For " +
            "more information about GLIMMPSE and related publications, please visit\n\n";
    private static final String TEXT_INTRO2 = "The automated validation tests shown below " +
            "compare power values produced by the JavaStatistics library to published results " +
            "and also to simulation.  Sources for published values include POWERLIB (Johnson \\emph{et al.} 2007) and " +
            "a SAS IML implementation of the methods described by Glueck and Muller (2003).\n\n" +
    		"Validation results are listed in Section 3 of the report.  Timing results show the calculation " +
    		"and simulation times for the overall experiment and the mean times per power calculation.  " +
    		"Summary statistics show the maximum absolute deviation between the power value calculated " +
    		"by the JavaStatistics library and the results obtained from SAS or via simulation.  The table in " +
    		"Section 3.3 shows the deviation values for each individual power comparison.  Deviations larger " +
    		"than $10^{-6}$ from SAS power values and $0.05$ for simulated power values are displayed in red.\n\n ";
    private static final String REFERENCES = "\\section*{References}\n\n" +
            "\\hangindent2em\n\\hangafter=1\n Glueck, D. H., \\& Muller, K. E. (2003). " +
            "Adjusting power for a baseline covariate in linear models. \\emph{Statistics " +
            "in Medicine}, \\emph{22}(16), 2535-2551.\n\n" +
            "\\hangindent2em\n\\hangafter=1\n Johnson, J. L., Muller, K. E., Slaughter, " +
            "J. C., Gurka, M. J., \\& Gribbin, M. J. (2009). POWERLIB: SAS/IML Software for Computing " +
            "Power in Multivariate Linear Models. \\emph{Journal of Statistical Software}, " +
            "\\emph{30}(5), 1-27.\n\n" +
            "\\hangindent2em\n\\hangafter=1\n Muller, K. E., \\& Stewart, P. W. (2006). " +
            "Linear model theory: univariate, multivariate, " +
            "and mixed models. Hoboken, New Jersey: John Wiley and Sons.";
    private static DecimalFormat Number = new DecimalFormat("#0.0000000");
    private static DecimalFormat ShortNumber = new DecimalFormat("#0.0000");
    private static DecimalFormat VeryShortNumber = new DecimalFormat("#0.00");
    private static DecimalFormat LongNumber = new DecimalFormat("#0.00000000");
    private static DecimalFormat ScientificNumber = new DecimalFormat("0.00E0");  

    // tolerance for comparison
    private double sasTolerance = 0.00001;
    private double simTolerance = 0.05;

    /**
     * Create a validation report builder using default
     * cutoffs for comparison against the calculated powers.
     */
    public ValidationReportBuilder() {}

    /**
     * Create a validation report builder with the specified
     * cutoffs for comparison against the calculated powers.
     * @param sasTolerance comparison threshold for SAS
     * @param simTolerance comparison threshold for simulation
     */
    public ValidationReportBuilder(double sasTolerance, 
            double simTolerance) {
        this.sasTolerance = sasTolerance;
        this.simTolerance = simTolerance;
    }


    /**
     * Write the validation report to LaTex.
     */
    public void createValidationReportAsLaTex(String filename,
            String title, String author, String studyDesignDescription,
            GLMMPowerParameters params,
            PowerChecker checker) 
                    throws FileNotFoundException, IOException {
        createValidationReportAsLaTex(filename, title, author, 
                studyDesignDescription, params, null, checker, null);
    }

    /**
     * Write the validation report to LaTex.
     */
    public void createValidationReportAsLaTex(String filename,
            String title, String author, String studyDesignDescription,
            GLMMPowerParameters params,
            PowerChecker checker, String image) 
                    throws FileNotFoundException, IOException {
        createValidationReportAsLaTex(filename, title, author, 
                studyDesignDescription, params, null, checker, image);
    }

    /**
     * Write the validation report to LaTex.
     */
    public void createValidationReportAsLaTex(String filename,
            String title, String author, String studyDesignDescription,
            GLMMPowerParameters params, String matrixAltString,
            PowerChecker checker) 
                    throws FileNotFoundException, IOException {
        createValidationReportAsLaTex(filename, title, author, 
                studyDesignDescription, params, matrixAltString, checker, null);
    }

    /**
     * Write the validation report to LaTex.
     */
    public void createValidationReportAsLaTex(String filename,
            String title, String author, String studyDesignDescription,
            GLMMPowerParameters params, String matrixAltString, 
            PowerChecker checker, String image) 
                    throws FileNotFoundException, IOException {
        if (filename != null) {
            FileWriter fWriter = null;
            BufferedWriter outFile = null;
            try {
                fWriter = new FileWriter(filename);
                outFile = new BufferedWriter(fWriter);
                addPreamble(outFile, title, author);
                addIntroduction(outFile, title, author);
                addStudyDesignInfo(outFile, studyDesignDescription, params, matrixAltString);
                addResults(outFile, checker, (params.getSigmaGaussianRandom() != null));
                if (image != null) {
                    outFile.write("\n\n\\includegraphics[width=4in]{" + image + "}\n\n");
                }
                outFile.write(REFERENCES);
                addClosing(outFile);
            } catch (Exception e) {
                throw new IOException(
                        "Failed to create validation report: " + e.getMessage());
            } finally {
                if (outFile != null) {
                    outFile.close();
                }
                if (fWriter != null) {
                    fWriter.close();
                }
            }
        }
    }    

    /**
     * Write the latex preamble for a pdflatex document
     * @param outFile output file stream
     * @throws IOException
     */
    private void addPreamble(BufferedWriter outFile,
            String title, String author) 
            throws IOException {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();
        
        outFile.write("\\documentclass{glimmpse-report}\n" +
                "\\JavaStatisticsVersion{" + LIBRARY_VERSION +"}\n" +
                "\\doctitle{" + title + "}\n" + 
                "\\docauthor{" + author + "}\n" + 
                "\\docdate{" + dateFormat.format(date) + "}\n");
        
//        outFile.write("\\documentclass[12pt,english]{article}\n" +
//                "\\renewcommand{\\familydefault}{\\sfdefault}" +
//                "\\usepackage{longtable,tabu}\n" +
//                "\\usepackage{savesym}\n" +
//                "\\usepackage{amsmath}\n" +
//                "\\savesymbol{iint}\n" +
//                "\\usepackage[T1]{fontenc}\n" +
//                "\\usepackage[latin9]{inputenc}\n" +
//                "\\usepackage{esint}\n" +
//                "\\usepackage{babel}\n" +
//                "\\usepackage{hyperref}\n" +
//                "\\usepackage{geometry}\n" +
//                "\\geometry{verbose,tmargin=1in,bmargin=1in,lmargin=0.5in,rmargin=0.5in}\n" +
//                "\\setlength{\\parskip}{\\medskipamount}\n" +
//                "\\setlength{\\parindent}{0pt}\n" +
//                "\\PassOptionsToPackage{normalem}{ulem}\n" +
//                "\\usepackage{ulem}\n" +
//                "\\usepackage{tabularx}\n" +
//                "\\usepackage{graphicx}\n" +
//                "\\usepackage{color}\n" +
//                "\\makeatletter\n\n" +
//                "\\begin{document}\n\n" +
//                "\\setcounter{MaxMatrixCols}{100}");
    }

    /**
     * Write the latex document closure.
     * @param outFile output file stream
     * @throws IOException
     */
    private void addClosing(BufferedWriter outFile) 
            throws IOException {
        outFile.write("\n\\end{document}\n");
    }

    /**
     * Add the introduction chapter.
     * @param outFile
     * @param title
     * @param author
     */
    private void addIntroduction(BufferedWriter outFile,
            String title, String author) 
                    throws IOException {
        // start document
        outFile.write("\\begin{document}\n");
        // add the introduction section shared across all validation reports
        outFile.write("\\section{"+SECTION_INTRO +"}\n");
        outFile.write(TEXT_INTRO1 + 
                " \n\n\\href{http://samplesizeshop.org}{http://samplesizeshop.org}.\n\n" +
                TEXT_INTRO2);
    }

    /**
     * Add the study design info
     * @param outFile
     * @param title
     * @param author
     */
    private void addStudyDesignInfo(BufferedWriter outFile,
            String studyDesignDescription, 
            GLMMPowerParameters params,
            String matrixAltString) 
                    throws IOException {

        outFile.write("\\section{"+SECTION_DESIGN +"}\n");
        outFile.write(studyDesignDescription + "\n");
        outFile.write("\\subsection{"+SUBSECTION_INPUTS +"}\n");
        // add the inputs
        addListInputs(outFile, params);
        addMatrixInputs(outFile, params, matrixAltString);   
    }

    /**
     * Add the results section.
     * @param outFile
     * @param checker
     * @throws IOException
     */
    private void addResults(BufferedWriter outFile,
            PowerChecker checker, boolean hasCovariate) 
                    throws IOException {
        outFile.write("\\section{"+SECTION_RESULTS +"}\n");
        outFile.write("A total of "+ checker.getResults().size() +
                " power values were computed for this experiment.\n\n");
        // add timing results
        outFile.write("\\subsection{"+SUBSECTION_TIMING +"}\n");
        addTimingResults(outFile, checker.getTiming(), checker.getResults().size());
        // add the summary stats
        outFile.write("\\subsection{"+SUBSECTION_SUMMARY +"}\n");
        addSummaryStatistics(outFile, checker);
        // add the results table
        outFile.write("\\subsection{"+SUBSECTION_RESULTS +"}\n");
        addResultsTable(outFile, checker, hasCovariate);
    }

    /**
     * Display the list inputs used in the study design.
     * @param params input parameters to the power calculation
     */
    private void addListInputs(BufferedWriter outFile,
            GLMMPowerParameters params) 
                    throws IOException {
        outFile.write("\\subsubsection{List Inputs}\n\n");
        writeDoubleList(outFile, "Type I error rates", params.getAlphaList());
        writeDoubleList(outFile, "Beta scale values", params.getBetaScaleList());
        writeDoubleList(outFile, "Sigma scale values", params.getSigmaScaleList());
        writeIntegerList(outFile, "Per group sample size values", params.getSampleSizeList());
        writeDoubleList(outFile, "Nominal power values", params.getPowerList());
        List<Test> testList = params.getTestList();
        if (testList != null && testList.size() > 0) {
            outFile.write("{\\bf Statistical tests}\n\n");
            boolean first = true;
            for(Test value: testList) {
                if (!first) {
                    outFile.write(", ");
                }
                outFile.write(testToString(value));
                if (first) {
                    first = false;
                }
            }
            outFile.write("\n\n");
        }
        List<PowerMethod> powerMethodList = params.getPowerMethodList();
        if (powerMethodList != null && powerMethodList.size() > 0) {
            outFile.write("{\\bf Power methods}\n\n");
            boolean first = true;
            for(PowerMethod value: powerMethodList) {
                if (!first) {
                    outFile.write(", ");
                }
                outFile.write(powerMethodToString(value));
                if (first) {
                    first = false;
                }
            }
            outFile.write("\n\n");
        }
        writeDoubleList(outFile, "Power quantile values:", params.getQuantileList());

    }  

    /**
     * Write out a list of integers
     * @param section
     * @param label
     * @param list
     */
    private void writeIntegerList(BufferedWriter outFile, 
            String label, List<Integer> list) 
                    throws IOException {
        if (list != null && list.size() > 0) {
            outFile.write("{\\bf " + label + "}\n\n");
            boolean first = true;
            for(Integer value: list) {
                if (!first) {
                    outFile.write(", ");
                }
                outFile.write(Integer.toString(value));
                if (first) {
                    first = false;
                }
            }
            outFile.write("\n\n");
        }
    }

    /**
     * Write out a list of integers
     * @param section
     * @param label
     * @param list
     */
    private void writeDoubleList(BufferedWriter outFile, 
            String label, List<Double> list) 
                    throws IOException {
        if (list != null && list.size() > 0) {
            outFile.write("{\\bf " + label + "}\n\n");
            boolean first = true;
            for(Double value: list) {
                if (!first) {
                    outFile.write(", ");
                }
                outFile.write(Number.format(value));
                if (first) {
                    first = false;
                }
            }
            outFile.write("\n\n");
        }
    }

    /**
     * Display the matrix inputs used in the study design.
     * @param outFile output LaTeX file
     * @param params input parameters to the power calculation
     * @param matrixAltString alternative matrix latex for large or
     * other difficult to typeset matrices 
     */
    private void addMatrixInputs(BufferedWriter outFile,
            GLMMPowerParameters params,
            String matrixAltString) 
                    throws IOException {
        outFile.write("\\subsubsection{Matrix Inputs}\n\n");

        if (matrixAltString != null) {
            outFile.write(matrixAltString + "\n\n");
        } else {
            RealMatrix xEssence = params.getDesignEssence();
            RealMatrix beta = params.getBeta().getCombinedMatrix();
            RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
            RealMatrix U = params.getWithinSubjectContrast();
            RealMatrix thetaNull = params.getTheta();
            RealMatrix sigmaE = params.getSigmaError();
            RealMatrix sigmaG = params.getSigmaGaussianRandom();
            RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
            RealMatrix sigmaY = params.getSigmaOutcome();

            if (xEssence != null) {
                writeMatrix(outFile, "\\text{Es}\\left(\\mathbf{X}\\right)", xEssence);
            }
            if (beta != null) {
                writeMatrix(outFile, "\\mathbf{B}", beta);
            }
            if (C != null) {
                writeMatrix(outFile, "\\mathbf{C}", C);
            }
            if (U != null) {
                writeMatrix(outFile, "\\mathbf{U}", U);
            }
            if (thetaNull != null) {
                writeMatrix(outFile, "\\mathbf{\\Theta}_0", thetaNull);
            }
            if (sigmaE != null) {
                writeMatrix(outFile, "\\mathbf{\\Sigma}_E", sigmaE);
            }
            if (sigmaY != null) {
                writeMatrix(outFile, "\\mathbf{\\Sigma}_Y", sigmaY);
            }
            if (sigmaG != null) {
                writeMatrix(outFile, "\\mathbf{\\Sigma}_g", sigmaG);
            }
            if (sigmaYG != null) {
                writeMatrix(outFile, "\\mathbf{\\Sigma}_Yg", sigmaYG);
            }
        }
    }   

    /**
     * Write a matrix in latex
     * @param section
     * @param name
     * @param matrix
     */
    private void writeMatrix(BufferedWriter outFile, String name,
            RealMatrix matrix) throws IOException {

        if (matrix.getColumnDimension() > 10) {
            outFile.write("\\tiny\n");
        }
        outFile.write("\\begin{eqnarray*}\n");
        // add name label
        outFile.write("\\underset{\\left("+ matrix.getRowDimension() + 
                        "\\times" + matrix.getColumnDimension() + "\\right)}{" + name + 
                        "} & = & \\begin{bmatrix}");
        for(int r = 0; r < matrix.getRowDimension(); r++) {
            boolean first = true;
            for(int c = 0; c < matrix.getColumnDimension(); c++) {
                if (!first) {
                    outFile.write(" & ");
                }
                if (matrix.getColumnDimension() > 10) {
                    outFile.write(VeryShortNumber.format(matrix.getEntry(r, c)));
                } else {
                    outFile.write(ShortNumber.format(matrix.getEntry(r, c)));
                }
                if (first) {
                    first = false;
                }
            }
            outFile.write("\\protect\\\\\n");
        }
        outFile.write("\\end{bmatrix}\n");
        outFile.write("\\end{eqnarray*}\n");
        if (matrix.getColumnDimension() > 10) {
            outFile.write("\\normalsize\n");
        }
    }

    /**
     * Add a section with timing results.
     * @param timer timer object from the power check
     */
    private void addTimingResults(BufferedWriter outFile,
            PowerChecker.Timer timer, int totalCases) 
                    throws IOException {

        double calcTimeTotalSeconds = ((double) timer.calculationMilliseconds / 1000.0);
        double calcTimeAvgSeconds = 
                ((double) timer.calculationMilliseconds / (double) totalCases) / 1000.0;
        double simTimeTotalSeconds = ((double) timer.simulationMilliseconds / 1000.0);
        double simTimeAvgSeconds = 
                ((double) timer.simulationMilliseconds / (double) totalCases) / 1000.0;

        outFile.write("\\begin{tabular}{|l|l|l|}\n");
        outFile.write("\\hline\n");
        outFile.write(" & Total Time (seconds) & Mean Time (seconds) \\\\ \n");
        outFile.write("\\hline\n");
        outFile.write("Calculation & "+ Number.format(calcTimeTotalSeconds) + " & " +
                ScientificNumber.format(calcTimeAvgSeconds) + "\\tabularnewline\n\\hline\n");
        outFile.write("Simulation & "+  Number.format(simTimeTotalSeconds) + " & " +
                ScientificNumber.format(simTimeAvgSeconds) + "\\tabularnewline\n");
        outFile.write("\\hline\n");
        outFile.write("\\end{tabular}\n");
    }

    /**
     * Add summary statistics for the validation experiment to the PDF
     * @param section document section object
     * @param checker power checker with summary information
     */
    private void addSummaryStatistics(BufferedWriter outFile,
            PowerChecker checker) 
                    throws IOException{
        outFile.write("\\begin{tabular}{|l|l|}\n");
        outFile.write("\\hline\n");
        outFile.write("Max deviation from SAS & "+ 
                LongNumber.format(checker.getMaxSasDeviation())+
                "\\tabularnewline\n\\hline\n\n");
        if (checker.getMaxSaslowerCIDeviation() > -1 &&
                checker.getMaxSasUpperCIDeviation() > -1) {
            outFile.write("Max deviation from lower CI limit & "+ 
                    LongNumber.format(checker.getMaxSaslowerCIDeviation())+
                    "\\tabularnewline\n\\hline\n\n");
            outFile.write("Max deviation from upper CI limit & "+ 
                    LongNumber.format(checker.getMaxSasUpperCIDeviation())+
                    "\\tabularnewline\n\\hline\n\n");
        }
        outFile.write("Max deviation from simulation & "+ 
                LongNumber.format(checker.getMaxSimDeviation())+
                "\\tabularnewline\n\\hline\n\n");
        outFile.write("\\end{tabular}\n");     
    }

    /**
     * Add a section with power comparison results.  The section also
     * includes a brief header describing how to interpret the results.
     * 
     * @param results result list from the power check
     */
    private void addResultsTable(BufferedWriter outFile,
            PowerChecker checker, boolean hasCovariate) throws IOException {
        List<PowerChecker.Result> checkerResults = checker.getResults();
        if (checkerResults != null) {
            boolean hasCI = (checkerResults.size() > 0 && 
                    checkerResults.get(0).calculatedPower.getConfidenceInterval() != null);
            if (hasCI) {
                outFile.write("\\scriptsize");
                if (hasCovariate) {
                    outFile.write("\\begin{longtabu}{|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|}\n");
                } else {
                    outFile.write("\\begin{longtabu}{|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|}\n");
                }
            } else {
                if (hasCovariate) {
                    outFile.write("\\scriptsize");
                    outFile.write("\\begin{longtabu}{|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|}\n");
                } else {
                    outFile.write("\\begin{longtabu}{|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|X[l]|}\n");
                }
            }
            outFile.write("\\hline\n");    
            outFile.write("{\\bf Power} & ");
            if (hasCI) {
                outFile.write("{\\bf CI} & ");
            }
            outFile.write("{\\bf SAS Power (deviation)} & ");
            if (hasCI) {
                outFile.write("{\\bf SAS CI (deviation)} & ");
            }
            outFile.write("{\\bf Sim Power (deviation)} & ");
            outFile.write("{\\bf Test} & ");
            outFile.write("{\\bf Sigma Scale} & ");
            outFile.write("{\\bf Beta Scale} & ");
            outFile.write("{\\bf Total N} & ");
            if (hasCovariate) {
                outFile.write("{\\bf Alpha} & ");
                outFile.write("{\\bf Method} & ");
                outFile.write("{\\bf Quantile}");
            } else {
                outFile.write("{\\bf Alpha} ");
            }
            outFile.write("\\\\ \\hline\n");


            for(Result result: checkerResults)
            {
                // add calculated power
                outFile.write(Number.format(result.calculatedPower.getActualPower()) + " & ");
                // if applicable, add confidence interval
                ConfidenceInterval ci = result.calculatedPower.getConfidenceInterval();
                if (ci != null) {
                    outFile.write("(" + Number.format(ci.getLowerLimit()) + ", " + 
                            Number.format(ci.getUpperLimit()) + ")" + " & ");
                }     
                // add SAS power and deviation 
                outFile.write(Number.format(result.sasPower) + " (");
                if (result.sasDeviation > sasTolerance) {
                    outFile.write("\\textcolor{red}{");
                }
                outFile.write(Number.format(result.sasDeviation));
                if (result.sasDeviation > sasTolerance) {
                    outFile.write("}) & ");
                } else {
                    outFile.write(") & ");
                }
                // if applicable, add SAS CI and deviation
                if (ci != null)
                {
                    outFile.write("(" + Number.format(result.sasCILower) + ", " + 
                            Number.format(result.sasCIUpper) + ") ");
                    outFile.write("\\{");
                    // lower ci deviation
                    if (result.sasCILowerDeviation > sasTolerance) {
                        outFile.write("\\textcolor{red}{");
                    }
                    outFile.write(Number.format(result.sasCILowerDeviation) + ", ");
                    if (result.sasCILowerDeviation > sasTolerance) {
                        outFile.write("}");
                    }
                    // upper ci deviation
                    if (result.sasCIUpperDeviation > sasTolerance) {
                        outFile.write("\\textcolor{red}{");
                    }
                    outFile.write(Number.format(result.sasCIUpperDeviation));
                    if (result.sasCIUpperDeviation > sasTolerance) {
                        outFile.write("\\}");
                    }
                    outFile.write("\\} & ");
                }
                // add simulation power
                outFile.write((Double.isNaN(result.simulatedPower) ? 
                        "N/A" :  Number.format(result.simulatedPower)) + 
                        " (");
                if (result.simulationDeviation > simTolerance) {
                    outFile.write("\\textcolor{red}{");
                }
                outFile.write((Double.isNaN(result.simulatedPower) ? 
                        "N/A" : Number.format(result.simulationDeviation)));
                if (result.simulationDeviation > simTolerance) {
                    outFile.write("}) & ");
                } else {
                    outFile.write(") & ");
                }

                outFile.write(testToString(result.calculatedPower.getTest()) + " & ");
                outFile.write(Number.format(result.calculatedPower.getSigmaScale()) + " & ");
                outFile.write(Number.format(result.calculatedPower.getBetaScale()) + " & ");
                outFile.write(Integer.toString(result.calculatedPower.getTotalSampleSize()) + " & ");
                if (hasCovariate) {
                    outFile.write(Number.format(result.calculatedPower.getAlpha()) + " & ");
                    outFile.write(powerMethodToString(result.calculatedPower.getPowerMethod()) + " & ");
                    outFile.write(Double.toString(result.calculatedPower.getQuantile()));
                } else {
                    outFile.write(Number.format(result.calculatedPower.getAlpha()));
                }
                outFile.write("\\\\ \\hline\n");
            }

            outFile.write("\\end{longtabu}\n");    
            outFile.write("\\normalsize\n");    
        }

    }

    /**
     * Write the results to an HTML file
     */
    public void createValidationReportAsHTML(PowerChecker checker,
            String title, String filename)
    {
        createValidationReportAsHTML(checker, title, filename, null);
    }

    /**
     * Write the results to an HTML file, optionally including a plot image
     * @param filename
     */
    public void createValidationReportAsHTML(PowerChecker checker,
            String title, String filename, String imageFilename)
    {
        List<PowerChecker.Result> checkerResults = checker.getResults();
        PowerChecker.Timer timer = checker.getTiming();
        // output the results
        StringBuffer buffer = new StringBuffer();

        buffer.append("<html><head></head><body><h1>" + title + "</h1>");
        buffer.append("<h3>Timing Results</h3>");
        buffer.append("<table border='1' cellpadding='5'>");
        buffer.append("<tr><td>Calculation time</td><td>" + ((double) timer.calculationMilliseconds / 1000)
                + "</td></tr>");
        buffer.append("<tr><td>Simulation time</td><td>" + ((double) timer.simulationMilliseconds / 1000)
                + "</td></tr></table><p></p>");

        buffer.append("<h3>Max Absolute Deviation Summary</h3>");
        buffer.append("<table border='1' cellpadding='5'>");
        buffer.append("<tr><td>Max deviation from SAS</td><td>"+
                LongNumber.format(checker.getMaxSasDeviation())+"</td></tr>");
        if (checker.getMaxSaslowerCIDeviation() >= 0 && checker.getMaxSasUpperCIDeviation() >= 0)
        {
            buffer.append("<tr><td>Max deviation from lower CI limit</td><td>"+
                    LongNumber.format(checker.getMaxSaslowerCIDeviation())+"</td></tr>");
            buffer.append("<tr><td>Max deviation from upper CI limit</td><td>"+
                    LongNumber.format(checker.getMaxSasUpperCIDeviation())+"</td></tr>");
        }
        buffer.append("<tr><td>Max deviation from simulation</td><td>"+
                LongNumber.format(checker.getMaxSimDeviation())+"</td></tr>");
        buffer.append("</table><p></p>");

        buffer.append("<h3>Full Results</h3>");
        buffer.append("<table border='1' cellpadding='5'><tr><th>Calc Power</th>");
        if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
        {
            ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
            buffer.append("<th>Confidence Interval (CI)</th><th>&alpha;-lower</th><th>&alpha;-upper</th>");
        }
        buffer.append("<th>SAS Power (deviation)</th>");
        if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
        {
            ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
            buffer.append("<th>SAS CI {deviation}</th>");
        }
        buffer.append("<th>Sim Power (deviation)</th>");
        buffer.append("<th>Test</th><th>Sigma Scale</th><th>Beta Scale</th><th>Total N</th>");
        buffer.append("<th>Alpha</th><th>Power Method</th><th>Quantile</th></tr>");

        for(Result result: checkerResults)
        {
            buffer.append("<tr><td>");
            buffer.append(Number.format(result.calculatedPower.getActualPower()));
            buffer.append("</td><td>");
            ConfidenceInterval ci = result.calculatedPower.getConfidenceInterval();
            if (ci != null)
            {
                buffer.append("(" + Number.format(ci.getLowerLimit()) + ", " + Number.format(ci.getUpperLimit()) + ")");
                buffer.append("</td><td>");
                buffer.append(ci.getAlphaLower());
                buffer.append("</td><td>");
                buffer.append(ci.getAlphaUpper());
                buffer.append("</td><td>");
            }
            buffer.append(Number.format(result.sasPower));
            if (result.sasDeviation > sasTolerance)
                buffer.append(" <font color='red'>(" + Number.format(result.sasDeviation) + ")</font></td><td>");
            else
                buffer.append(" (" + Number.format(result.sasDeviation) + ")</td><td>");
            if (ci != null)
            {
                buffer.append("(" + Number.format(result.sasCILower) + ", " + Number.format(result.sasCIUpper) + ")");
                buffer.append("{"); 
                if (result.sasCILowerDeviation > sasTolerance)
                    buffer.append("<font color='red'>" + Number.format(result.sasCILowerDeviation) + "</font>");
                else
                    buffer.append(Number.format(result.sasCILowerDeviation));
                buffer.append(", ");
                if (result.sasCIUpperDeviation > sasTolerance)
                    buffer.append("<font color='red'>" + Number.format(result.sasCIUpperDeviation) + "</font>");
                else
                    buffer.append(Number.format(result.sasCIUpperDeviation));
                buffer.append("}</td><td>");
            }
            buffer.append(Number.format(result.simulatedPower));
            if (result.simulationDeviation > simTolerance)
                buffer.append(" <font color='red'>(" + Number.format(result.simulationDeviation) + ")</font></td><td>");
            else
                buffer.append(" (" + Number.format(result.simulationDeviation) + ")</td><td>");

            buffer.append(result.calculatedPower.getTest() + "</td><td>" + 
                    Number.format(result.calculatedPower.getSigmaScale()) + "</td><td>" + 
                    Number.format(result.calculatedPower.getBetaScale()) + "</td><td>" + 
                    result.calculatedPower.getTotalSampleSize() + "</td><td>" + 
                    Number.format(result.calculatedPower.getAlpha()) + "</td><td>" + 
                    powerMethodToString(result.calculatedPower.getPowerMethod()) + "</td><td>" + 
                    result.calculatedPower.getQuantile() + "</td></tr>");
        }

        buffer.append("</table><p>");

        if (imageFilename != null)
        {
            buffer.append("<p><img src='" + imageFilename + "' /></p>");
        }

        buffer.append("</body></html>");

        FileWriter writer = null;
        try
        {
            writer = new FileWriter(filename);
            writer.write(buffer.toString());
        }
        catch (Exception e)
        {
            // TODO:
        }
        finally
        {
            if (writer != null) try {writer.close(); } catch (Exception e) {};
        }

    }

    /**
     * Write the current result set to stdout.
     */
    public void createValidationReportAsStdout(PowerChecker checker, 
            String title, boolean verbose)
    {
        List<PowerChecker.Result> checkerResults = checker.getResults();
        PowerChecker.Timer timer = checker.getTiming();
        if (verbose)
        {
            // output the results
            System.out.println("Calc Power (lower, upper)\tSAS Power (deviation)\t" +
                    "Sim Power (deviation)\tTest\tSigmaScale\tBetaScale\tTotal N\tAlpha\tPowerMethod\tQuantile");

            for(Result result: checkerResults)
            {
                System.out.println(Number.format(result.calculatedPower.getActualPower()) + "(" + 
                        (result.calculatedPower.getConfidenceInterval() != null ? 
                                Number.format(result.calculatedPower.getConfidenceInterval().getLowerLimit()) : "n/a") + ", " + 
                                (result.calculatedPower.getConfidenceInterval() != null ? 
                                        Number.format(result.calculatedPower.getConfidenceInterval().getUpperLimit()) : "n/a") + ")\t" +  
                                        Number.format(result.sasPower) + " (" +  Number.format(result.sasDeviation) + ")\t" +
                                        Number.format(result.simulatedPower) + " (" + Number.format(result.simulationDeviation) + ")\t" +
                                        result.calculatedPower.getTest() + "\t" + 
                                        Number.format(result.calculatedPower.getSigmaScale()) + "\t" + 
                                        Number.format(result.calculatedPower.getBetaScale()) + "\t" + 
                                        result.calculatedPower.getTotalSampleSize() + "\t" + 
                                        Number.format(result.calculatedPower.getAlpha()) + "\t" + 
                                        result.calculatedPower.getPowerMethod() + "\t" + 
                                        result.calculatedPower.getQuantile() + "\t");
            }
        }
        // output summary information
        System.out.println("===========================================");
        System.out.println("Summary statistics: " + title);
        System.out.println("===========================================");
        System.out.println("Total Calculation Time: " + ((double) timer.calculationMilliseconds / 1000));
        System.out.println("Mean Calculation Time: " + (((double) timer.calculationMilliseconds / 1000)/checkerResults.size()));
        System.out.println("Total Simulation Time: " + ((double) timer.simulationMilliseconds / 1000));
        System.out.println("Mean Simulation Time: " + (((double) timer.simulationMilliseconds / 1000)/checkerResults.size()));
        System.out.println("Max Deviation from Published: " + LongNumber.format(checker.getMaxSasDeviation()));
        System.out.println("Max Deviation from Simulation: " + LongNumber.format(checker.getMaxSimDeviation()));
        if (checker.getMaxSaslowerCIDeviation() >= 0 || checker.getMaxSasUpperCIDeviation() >=0)
        {
            System.out.println("Max Deviation from SAS, Lower Confidence Limit: " + 
                    LongNumber.format(checker.getMaxSaslowerCIDeviation()));
            System.out.println("Max Deviation from SAS, Upper Confidence Limit: " + 
                    LongNumber.format(checker.getMaxSasUpperCIDeviation()));
        }
    }

    /**
     * Pretty display of power method.
     * @param method power computation method
     * @return display string
     */
    private String powerMethodToString(GLMMPowerParameters.PowerMethod method)
    {
        String value = null;
        switch (method)
        {
        case CONDITIONAL_POWER:
            value = "cond";
            break;
        case UNCONDITIONAL_POWER:
            value = "uncond";
            break;
        case QUANTILE_POWER:
            value = "quantile";
            break;          
        }
        return value;
    }

    /**
     * Pretty display of statistical test.
     * @param method power computation method
     * @return display string
     */
    private String testToString(GLMMTestFactory.Test test)
    {
        String value = null;
        switch (test)
        {
        case HOTELLING_LAWLEY_TRACE:
            value = "HLT";
            break;
        case WILKS_LAMBDA:
            value = "WL";
            break;
        case PILLAI_BARTLETT_TRACE:
            value = "PBT";
            break;          
        case UNIREP:
            value = "UNIREP";
            break;   
        case UNIREP_GEISSER_GREENHOUSE:
            value = "UNIREP-GG";
            break;   
        case UNIREP_HUYNH_FELDT:
            value = "UNIREP-HF";
            break;   
        case UNIREP_BOX:
            value = "UNIREP-BOX";
            break;   
        }
        return value;
    }

}
