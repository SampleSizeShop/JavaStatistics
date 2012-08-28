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

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math.linear.RealMatrix;

import com.itextpdf.text.Anchor;
import com.itextpdf.text.BaseColor;
import com.itextpdf.text.Chapter;
import com.itextpdf.text.Chunk;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Font;
import com.itextpdf.text.Paragraph;
import com.itextpdf.text.Section;
import com.itextpdf.text.pdf.PdfPTable;
import com.itextpdf.text.pdf.PdfWriter;

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
    private static int PARAGRAPH_SPACING = 12;
    private static Font TITLE_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 18,
            Font.BOLD);
    private static Font SECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 16,
            Font.BOLD);
    private static Font SUBSECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 14,
            Font.BOLD);
    private static Font URL_FONT= new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.NORMAL, BaseColor.BLUE);
    private static Font UNDERLINE_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.UNDERLINE);
    private static Font RED_FONT= new Font(Font.FontFamily.TIMES_ROMAN, 8,
            Font.NORMAL, BaseColor.RED);
    private static Font BOLD_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 8,
            Font.BOLD);
    private static Font TABLE_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 8,
            Font.NORMAL);
    
    // section headers and static text
    private static final String TITLE_PREFIX = "Test Case Validation Report: ";
    private static final String SAMPLESIZESHOP_URL = "http://samplesizeshop.com";
    private static final String SECTION_INTRO = "Introduction";
    private static final String SECTION_DESIGN = "Study Design";
    private static final String SECTION_RESULTS = "Validation Results";   
    private static final String SUBSECTION_INPUTS = "Inputs to the Power Calculation";
    private static final String SUBSECTION_TIMING = "Timing";
    private static final String SUBSECTION_SUMMARY = "Summary Statistics";
    private static final String SUBSECTION_RESULTS = "Full Validation Results";
    private static final String TEXT_INTRO1 = "The following report contains validation results " +
            "for the JavaStatistics library, a component of the GLIMMPSE software system.  For " +
            "more information about GLIMMPSE and related publications, please visit ";
    private static final String TEXT_INTRO2 = "The automated validation tests shown below " +
            "compare power values produced by the JavaStatistics library against published results" +
            "and simulation.  Sources for published values include POWERLIB (Johnson 2007) and " +
            "a SAS IML implementations of the methods described by Glueck & Muller (2003).";

    private static DecimalFormat Number = new DecimalFormat("#0.0000000");
    private static DecimalFormat ShortNumber = new DecimalFormat("#0.0000");
    private static DecimalFormat LongNumber = new DecimalFormat("#0.00000000");

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
     * Open a new pdf document for writing.
     * 
     * @param filename path to pdf file
     * @return pdf document object
     * @throws DocumentException
     * @throws FileNotFoundException
     */
    public Document openDocument(String filename) 
            throws DocumentException, FileNotFoundException {
        Document document = new Document();
        PdfWriter.getInstance(document, new FileOutputStream(filename));
        document.open();
        return document;
    }

    /**
     * Close the pdf document.
     * @param document the pdf document object
     */
    public void closeDocument(Document document) {
        if (document != null) {
            document.close();
        }
    }

    /**
     * Write the validation report to LaTex.
     * @param document the pdf document
     * @param timer timer object from the power check
     */
    public void createValidationReport(String filename,
            String title, String author, String studyDesignDescription,
            GLMMPowerParameters params,
            PowerChecker checker) 
                    throws DocumentException {
        Document document = null;
        try {
            document = openDocument(filename);
            addIntroduction(document, title, author);
            addStudyDesignInfo(document, studyDesignDescription, params);
            addResults(document, checker);
        } catch (Exception e) {
            throw new DocumentException(
                    "Failed to create validation report: " + e.getMessage());
        } finally {
            if (document != null) {
                document.close();
            }
        }
    }    

    /**
     * Add document meta data and add the introduction chapter.
     * @param document
     * @param title
     * @param author
     */
    private void addIntroduction(Document document,
            String title, String author) 
                    throws DocumentException {
        // add meta data
        document.addTitle(TITLE_PREFIX + title);
        document.addSubject(
                "Validation results for the JavaStatistics library in " +
                "the GLIMMPSE software system");
        document.addKeywords("GLIMMPSE, Java, statistics");
        if (author != null) {
            document.addAuthor(author);
            document.addCreator(author);
        }
        // add title
        Paragraph preface = createParagraph(null, null);
        preface.add(new Paragraph(TITLE_PREFIX, TITLE_FONT));
        preface.add(new Paragraph(title, TITLE_FONT));

        // add the introduction section shared across all validation reports
        Chapter chapter = createChapter(SECTION_INTRO, 1);
        Paragraph introP1 = createParagraph(null, null);
        Anchor url = new Anchor(SAMPLESIZESHOP_URL + ".", URL_FONT);
        url.setReference(SAMPLESIZESHOP_URL);
        introP1.add(TEXT_INTRO1);
        introP1.add(url);
        chapter.add(introP1);
        chapter.add(createParagraph(TEXT_INTRO2, null));

        // add to the page
        document.add(preface);
        document.add(chapter);
    }


    private void addStudyDesignInfo(Document document,
            String studyDesignDescription, 
            GLMMPowerParameters params) 
                    throws DocumentException {
        Chapter chapter = createChapter(SECTION_DESIGN, 2);
        // add the study design description
        chapter.add(createParagraph(studyDesignDescription, null));
        // add the inputs
        Section inputsSection = createSection(chapter, SUBSECTION_INPUTS);
        addListInputs(inputsSection, params);
        addMatrixInputs(inputsSection, params);   

        document.add(chapter);
    }

    private void addResults(Document document,
            PowerChecker checker) 
                    throws DocumentException {
        Chapter chapter = createChapter(SECTION_RESULTS, 3);
        // add the timing results
        Section timingSection = createSection(chapter, SUBSECTION_TIMING);
        addTimingResults(timingSection, checker.getTiming());
        // add the summary stats
        Section summarySection = createSection(chapter, SUBSECTION_SUMMARY);
        addSummaryStatistics(summarySection, checker);
        // add the results table
        Section resultsSection = createSection(chapter, SUBSECTION_RESULTS);
        addResultsTable(resultsSection, checker);

        document.add(chapter);
    }

    /**
     * Display the list inputs used in the study design.
     * @param params input parameters to the power calculation
     */
    private void addListInputs(Section section,
            GLMMPowerParameters params) {

        writeDoubleList(section, "Type I error rates", params.getAlphaList());
        writeDoubleList(section, "Beta scale values", params.getBetaScaleList());
        writeDoubleList(section, "Sigma scale values", params.getSigmaScaleList());
        writeIntegerList(section, "Per group sample size values", params.getSampleSizeList());
        writeDoubleList(section, "Nominal power values", params.getPowerList());
        List<Test> testList = params.getTestList();
        if (testList != null) {
            Paragraph para = createParagraph(null,null);
            para.add(new Chunk("Statistical tests", UNDERLINE_FONT));
            para.add(Chunk.NEWLINE);
            boolean first = true;
            for(Test value: testList) {
                if (!first) {
                    para.add(", ");
                }
                para.add(value.toString());
                if (first) {
                    first = false;
                }
            }
            section.add(para);
        }
        List<PowerMethod> powerMethodList = params.getPowerMethodList();
        if (powerMethodList != null) {
            Paragraph para = createParagraph(null,null);
            para.add(new Chunk("Power methods", UNDERLINE_FONT));
            para.add(Chunk.NEWLINE);
            boolean first = true;
            for(PowerMethod value: powerMethodList) {
                if (!first) {
                    para.add(", ");
                }
                para.add(powerMethodToString(value));
                if (first) {
                    first = false;
                }
            }
            section.add(para);
        }
        writeDoubleList(section, "Power quantile values:", params.getQuantileList());

    }  

    /**
     * Write out a list of integers
     * @param section
     * @param label
     * @param list
     */
    private void writeIntegerList(Section section, String label, List<Integer> list) {
        if (list != null && list.size() > 0) {
            Paragraph para = createParagraph(null,null);
            para.add(new Chunk(label, UNDERLINE_FONT));
            para.add(Chunk.NEWLINE);
            boolean first = true;
            for(Integer value: list) {
                if (!first) {
                    para.add(", ");
                }
                para.add(Integer.toString(value));
                if (first) {
                    first = false;
                }
            }
            section.add(para);
        }
    }

    /**
     * Write out a list of integers
     * @param section
     * @param label
     * @param list
     */
    private void writeDoubleList(Section section, String label, List<Double> list) {
        if (list != null & list.size() > 0) {
            Paragraph para = createParagraph(null,null);
            para.add(new Chunk(label, UNDERLINE_FONT));
            para.add(Chunk.NEWLINE);
            boolean first = true;
            for(Double value: list) {
                if (!first) {
                    para.add(", ");
                }
                para.add(ShortNumber.format(value));
                if (first) {
                    first = false;
                }
            }
            section.add(para);
        }
    }

    /**
     * Display the matrix inputs used in the study design.
     * @param params input parameters to the power calculation
     */
    private void addMatrixInputs(Section section,
            GLMMPowerParameters params) {
        RealMatrix xEssence = params.getDesignEssence();
        RealMatrix beta = params.getBeta().getCombinedMatrix();
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getWithinSubjectContrast();
        RealMatrix thetaNull = params.getTheta();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();

        PdfPTable table = new PdfPTable(2);

        if (xEssence != null) {
            writeMatrix(section, "Design essence", xEssence);
        }
        if (beta != null) {
            writeMatrix(section, "Beta", beta);
        }
        if (C != null) {
            writeMatrix(section, "Between participant contrast", C);
        }
        if (U != null) {
            writeMatrix(section, "Within participant contrast", U);
        }
        if (thetaNull != null) {
            writeMatrix(section, "Null hypothesis", thetaNull);
        }
        if (sigmaE != null) {
            writeMatrix(section, "Covariance of errors", sigmaE);
        }
        if (sigmaY != null) {
            writeMatrix(section, "Covariance of outcomes", sigmaY);
        }
        if (sigmaG != null) {
            writeMatrix(section, "Covariance of Gaussian covariate", sigmaG);
        }
        if (sigmaYG != null) {
            writeMatrix(section, "Covariance of Gaussian covariate and outcomes", sigmaYG);
        }
    }   

    private void writeMatrix(Section section, String name,
            RealMatrix matrix) {
        // add name label
        section.add(createParagraph(name, UNDERLINE_FONT));
        // add matrix
        PdfPTable matrixTable = new PdfPTable(matrix.getColumnDimension());
        for(int r = 0; r < matrix.getRowDimension(); r++) {
            for(int c = 0; c < matrix.getColumnDimension(); c++) {
                matrixTable.addCell(ShortNumber.format(matrix.getEntry(r, c)));
            }
        }
        section.add(matrixTable);
    }

    /**
     * Create a new chapter with the specified label and number
     * @param label the chapter name/label
     * @param number the chapter number
     * @return Chapter object
     */
    private Chapter createChapter(String label, int number) {
        Paragraph paragraph = new Paragraph(label, 
                SECTION_FONT);
        paragraph.setSpacingAfter(PARAGRAPH_SPACING);
        Chapter chapter = new Chapter(paragraph, number);
        chapter.setTriggerNewPage(false);
        return chapter;
    }

    /**
     * Create a new section for the specified chapter
     * @param label the chapter name/label
     * @param number the chapter number
     * @return Chapter object
     */
    private Section createSection(Chapter chapter, String label) {
        Paragraph paragraph = new Paragraph(label, SUBSECTION_FONT);
        paragraph.setSpacingAfter(PARAGRAPH_SPACING);
        Section section = chapter.addSection(paragraph);
        return section;
    }

    /**
     * Create an empty paragraph with appropriate spacing.
     * @param contents optional contents for the paragraph
     * @param font optional font for the paragraph
     * @return Paragraph object
     */
    private Paragraph createParagraph(String contents, Font font) {
        Paragraph paragraph = new Paragraph();
        if (font != null) {
            paragraph.setFont(font);
        }
        paragraph.setSpacingAfter(PARAGRAPH_SPACING);
        if (contents != null) {
            paragraph.add(contents);
        }
        return paragraph;
    }

    /**
     * Add a section with timing results.
     * @param timer timer object from the power check
     */
    private void addTimingResults(Section section,
            PowerChecker.Timer timer) {
        PdfPTable table = new PdfPTable(2);
        // t.setBorderColor(BaseColor.GRAY);
        // t.setPadding(4);
        // t.setSpacing(4);
        // t.setBorderWidth(1);
        table.addCell("Calculation time (seconds)");
        table.addCell(Double.toString((double) timer.calculationMilliseconds / 1000.0));
        table.addCell("Simulation time (seconds)");
        table.addCell(Double.toString((double) timer.simulationMilliseconds / 1000.0));

        section.add(table);
    }

    /**
     * Add summary statistics for the validation experiment to the PDF
     * @param section document section object
     * @param checker power checker with summary information
     */
    private void addSummaryStatistics(Section section,
            PowerChecker checker) {

        PdfPTable table = new PdfPTable(2);
        table.addCell("Max deviation from SAS");
        table.addCell(LongNumber.format(checker.getMaxSasDeviation()));
        table.addCell("Max deviation from lower CI limit");
        table.addCell(LongNumber.format(checker.getMaxSaslowerCIDeviation()));
        table.addCell("Max deviation from upper CI limit");
        table.addCell(LongNumber.format(checker.getMaxSasUpperCIDeviation()));
        table.addCell("Max deviation from simulation");
        table.addCell(LongNumber.format(checker.getMaxSimDeviation()));
        section.add(table);        
    }

    /**
     * Add a section with power comparison results.  The section also
     * includes a brief header describing how to interpret the results.
     * 
     * @param results result list from the power check
     */
    private void addResultsTable(Section section,
            PowerChecker checker) {
        List<PowerChecker.Result> checkerResults = checker.getResults();
        if (checkerResults != null) {
            PdfPTable table = null;
            boolean hasCI = (checkerResults.size() > 0 && 
                    checkerResults.get(0).calculatedPower.getConfidenceInterval() != null);
            if (hasCI) {
                table = new PdfPTable(12);
            } else {
                table = new PdfPTable(10);
            }
            table.setWidthPercentage(100);
            table.addCell(createParagraph("Power", BOLD_FONT));
            if (hasCI) {
                table.addCell(createParagraph("CI", BOLD_FONT));
            }
            table.addCell(createParagraph("SAS Power (deviation)", BOLD_FONT));
            if (hasCI) {
                table.addCell(createParagraph("SAS CI (deviation)", BOLD_FONT));
            }
            table.addCell(createParagraph("Sim Power (deviation)", BOLD_FONT));
            table.addCell(createParagraph("Test", BOLD_FONT));
            table.addCell(createParagraph("Sigma Scale", BOLD_FONT));
            table.addCell(createParagraph("Beta Scale", BOLD_FONT));
            table.addCell(createParagraph("Total N", BOLD_FONT));
            table.addCell(createParagraph("Alpha", BOLD_FONT));
            table.addCell(createParagraph("Method", BOLD_FONT));
            table.addCell(createParagraph("Quantile", BOLD_FONT));

            for(Result result: checkerResults)
            {
                // add calculated power
                table.addCell(new Paragraph(
                        Number.format(result.calculatedPower.getActualPower()),
                        TABLE_FONT));
                // if applicable, add confidence interval
                ConfidenceInterval ci = result.calculatedPower.getConfidenceInterval();
                if (ci != null) {
                    table.addCell(new Paragraph("(" + 
                            Number.format(ci.getLowerLimit()) + ", " + 
                            Number.format(ci.getUpperLimit()) + ")", TABLE_FONT));
                }               
                // add the SAS power and deviation
                Paragraph sasCell = new Paragraph();
                sasCell.setFont(TABLE_FONT);
                sasCell.add(Number.format(result.sasPower));
                sasCell.add(Chunk.NEWLINE);
                Chunk sasDevChunk = new Chunk("(" + Number.format(result.sasDeviation) + ")");
                if (result.sasDeviation > sasTolerance) {
                    sasDevChunk.setFont(RED_FONT);
                }
                sasCell.add(sasDevChunk);
                table.addCell(sasCell);
                // if applicable, add SAS CI and deviation
                if (ci != null)
                {
                    Paragraph sasCICell = new Paragraph();
                    sasCICell.setFont(TABLE_FONT);
                    sasCICell.add(
                            "(" + Number.format(result.sasCILower) + ", " + 
                                    Number.format(result.sasCIUpper) + ")");
                    sasCICell.add(Chunk.NEWLINE);
                    sasCICell.add("{");
                    // lower ci deviation
                    Chunk sasCILowerDevChunk = new Chunk(
                            Number.format(result.sasCILowerDeviation));
                    if (result.sasCILowerDeviation > sasTolerance) {
                        sasCILowerDevChunk.setFont(RED_FONT);
                    }
                    sasCICell.add(sasCILowerDevChunk);
                    // upper ci deviation
                    Chunk sasCIUpperDevChunk = new Chunk(
                            Number.format(result.sasCIUpperDeviation));
                    if (result.sasCIUpperDeviation > sasTolerance) {
                        sasCIUpperDevChunk.setFont(RED_FONT);
                    }
                    sasCICell.add(sasCIUpperDevChunk);
                    sasCICell.add("}");
                }

                // add simulation power
                Paragraph simCell = new Paragraph();
                simCell.setFont(TABLE_FONT);
                simCell.add(Number.format(result.simulatedPower));
                simCell.add(Chunk.NEWLINE);
                Chunk simDevChunk = new Chunk("(" + Number.format(result.simulationDeviation) + ")");
                if (result.simulationDeviation > simTolerance) {
                    simDevChunk.setFont(RED_FONT);
                }
                simCell.add(simDevChunk);
                table.addCell(simCell);

                // add test
                table.addCell(new Paragraph(
                        result.calculatedPower.getTest().toString(), TABLE_FONT));
                table.addCell(new Paragraph(
                        Number.format(result.calculatedPower.getSigmaScale()), TABLE_FONT));
                table.addCell(new Paragraph(
                        Number.format(result.calculatedPower.getBetaScale()), TABLE_FONT));
                table.addCell(new Paragraph(
                        Integer.toString(result.calculatedPower.getTotalSampleSize()), TABLE_FONT));
                table.addCell(new Paragraph(
                        Number.format(result.calculatedPower.getAlpha()), TABLE_FONT));
                table.addCell(new Paragraph(
                        powerMethodToString(result.calculatedPower.getPowerMethod()), TABLE_FONT));
                table.addCell(new Paragraph(
                        Double.toString(result.calculatedPower.getQuantile()), TABLE_FONT));
            }
            section.add(table);
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

}
