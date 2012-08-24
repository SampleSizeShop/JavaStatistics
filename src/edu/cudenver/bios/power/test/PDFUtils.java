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
import com.itextpdf.text.Phrase;
import com.itextpdf.text.Section;
import com.itextpdf.text.pdf.PdfPCell;
import com.itextpdf.text.pdf.PdfPTable;
import com.itextpdf.text.pdf.PdfWriter;

import edu.cudenver.bios.power.glmm.GLMMTestFactory.Test;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters.PowerMethod;
import edu.cudenver.bios.power.test.PowerChecker.Result;
import edu.cudenver.bios.utils.ConfidenceInterval;

/**
 * Utility routines for creating PDF output for validation
 * results.  Depends on the iText library.
 * @author Sarah Kreidler
 *
 */
public class PDFUtils {
    private static int PARAGRAPH_SPACING = 12;
    private static Font TITLE_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 18,
            Font.BOLD);
    private static Font SECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 16,
            Font.BOLD);
    private static Font SUBSECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 14,
            Font.BOLD);
    private static Font URL_FONT= new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.NORMAL, BaseColor.BLUE);
    private static Font BOLD_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.BOLD);
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
    private static DecimalFormat LongNumber = new DecimalFormat("#0.00000000");
    
    /**
     * Open a new pdf document for writing.
     * 
     * @param filename path to pdf file
     * @return pdf document object
     * @throws DocumentException
     * @throws FileNotFoundException
     */
    public static Document openDocument(String filename) 
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
    public static void closeDocument(Document document) {
        if (document != null) {
            document.close();
        }
    }

    /**
     * Add the contents of the validation report to the
     * pdf.
     * @param document the pdf document
     * @param timer timer object from the power check
     */
    public static void createValidationReport(String filename,
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
    private static void addIntroduction(Document document,
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
    
    
    private static void addStudyDesignInfo(Document document,
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
    
    private static void addResults(Document document,
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
     * Display the matrix inputs used in the study design.
     * @param params input parameters to the power calculation
     */
    private static void addListInputs(Section section,
            GLMMPowerParameters params) {
        List<Double> alphaList = params.getAlphaList();
        List<Double> betaScaleList = params.getBetaScaleList();
        List<Double> sigmaScaleList = params.getSigmaScaleList();
        List<Integer> samplesSizeList = params.getSampleSizeList();
        List<Double> powerList = params.getPowerList();
        List<Test> testList = params.getTestList();
        List<PowerMethod> powerMethod = params.getPowerMethodList();
        List<Double> quantileList = params.getQuantileList();
        
        
        if (alphaList != null) {
            
        }
        if (betaScaleList != null) {
            
        }
        if (sigmaScaleList != null) {
            
        }
        if (samplesSizeList != null) {
            
        }
        if (powerList != null) {
            
        }
        if (testList != null) {
            
        }
        if (powerMethod != null) {
            
        }
        if (quantileList != null) {
            
        }
    }  
    
    /**
     * Display the matrix inputs used in the study design.
     * @param params input parameters to the power calculation
     */
    private static void addMatrixInputs(Section section,
            GLMMPowerParameters params) {
        RealMatrix xEssence = params.getDesignEssence();
        RealMatrix beta = params.getBeta().getCombinedMatrix();
        RealMatrix C = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix U = params.getBetweenSubjectContrast().getCombinedMatrix();
        RealMatrix thetaNull = params.getTheta();
        RealMatrix sigmaE = params.getSigmaError();
        RealMatrix sigmaG = params.getSigmaGaussianRandom();
        RealMatrix sigmaYG = params.getSigmaOutcomeGaussianRandom();
        RealMatrix sigmaY = params.getSigmaOutcome();
        
        PdfPTable table = new PdfPTable(2);

        if (xEssence != null) {
            section.add(createMatrixTable("Design essence", xEssence));
        }
        if (beta != null) {
            section.add(createMatrixTable("Beta", beta));
        }
        if (C != null) {
            section.add(createMatrixTable("Between participant contrast", C));
        }
        if (U != null) {
            section.add(createMatrixTable("Within participant contrast", U));
        }
        if (thetaNull != null) {
            section.add(createMatrixTable("Null hypothesis", thetaNull));
        }
        if (sigmaE != null) {
            section.add(createMatrixTable("Covariance of errors", sigmaE));
        }
        if (sigmaY != null) {
            section.add(createMatrixTable("Covariance of outcomes", sigmaY));
        }
        if (sigmaG != null) {
            section.add(createMatrixTable("Covariance of Gaussian covariate", sigmaG));
        }
        if (sigmaYG != null) {
            section.add(createMatrixTable("Covariance of Gaussian covariate and outcomes", sigmaYG));
        }
                
        section.add(table);

    }   
    
    private static PdfPTable createMatrixTable(String name,
            RealMatrix matrix) {
        PdfPTable table = new PdfPTable(2);
        PdfPCell cell = new PdfPCell();
        // add name 
        cell.setBorder(PdfPCell.NO_BORDER);
        cell.addElement(new Paragraph(name + " = "));
        table.addCell(cell);
        // add matrix
        PdfPTable matrixTable = new PdfPTable(matrix.getColumnDimension());
        for(int r = 0; r < matrix.getRowDimension(); r++) {
            for(int c = 0; c < matrix.getColumnDimension(); c++) {
                matrixTable.addCell(Number.format(matrix.getEntry(r, c)));
            }
        }
        PdfPCell matrixCell = new PdfPCell();
        matrixCell.setBorder(PdfPCell.NO_BORDER);
        matrixCell.addElement(matrixTable);
        table.addCell(matrixCell);
        return table;
    }
    
    /**
     * Create a new chapter with the specified label and number
     * @param label the chapter name/label
     * @param number the chapter number
     * @return Chapter object
     */
    private static Chapter createChapter(String label, int number) {
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
    private static Section createSection(Chapter chapter, String label) {
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
    private static Paragraph createParagraph(String contents, Font font) {
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
    private static void addTimingResults(Section section,
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
    private static void addSummaryStatistics(Section section,
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
    private static void addResultsTable(Section section,
            PowerChecker checker) {
        List<PowerChecker.Result> checkerResults = checker.getResults();
        if (checkerResults != null) {
            PdfPTable table = new PdfPTable(3);
//        buffer.append("<h3>Full Results</h3>");
//        buffer.append("<table border='1' cellpadding='5'><tr><th>Calc Power</th>");
//        if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
//        {
//            ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
//            buffer.append("<th>Confidence Interval (CI)</th><th>&alpha;-lower</th><th>&alpha;-upper</th>");
//        }
//        buffer.append("<th>SAS Power (deviation)</th>");
//        if (checkerResults.size() > 0 && checkerResults.get(0).calculatedPower.getConfidenceInterval() != null)
//        {
//            ConfidenceInterval ci = checkerResults.get(0).calculatedPower.getConfidenceInterval();
//            buffer.append("<th>SAS CI {deviation}</th>");
//        }
//        buffer.append("<th>Sim Power (deviation)</th>");
//        buffer.append("<th>Test</th><th>Sigma Scale</th><th>Beta Scale</th><th>Total N</th>");
//        buffer.append("<th>Alpha</th><th>Power Method</th><th>Quantile</th></tr>");
//
//        for(Result result: checkerResults)
//        {
//            buffer.append("<tr><td>");
//            buffer.append(Number.format(result.calculatedPower.getActualPower()));
//            buffer.append("</td><td>");
//            ConfidenceInterval ci = result.calculatedPower.getConfidenceInterval();
//            if (ci != null)
//            {
//                buffer.append("(" + Number.format(ci.getLowerLimit()) + ", " + Number.format(ci.getUpperLimit()) + ")");
//                buffer.append("</td><td>");
//                buffer.append(ci.getAlphaLower());
//                buffer.append("</td><td>");
//                buffer.append(ci.getAlphaUpper());
//                buffer.append("</td><td>");
//            }
//            buffer.append(Number.format(result.sasPower));
//            if (result.sasDeviation > sasTolerance)
//                buffer.append(" <font color='red'>(" + Number.format(result.sasDeviation) + ")</font></td><td>");
//            else
//                buffer.append(" (" + Number.format(result.sasDeviation) + ")</td><td>");
//            if (ci != null)
//            {
//                buffer.append("(" + Number.format(result.sasCILower) + ", " + Number.format(result.sasCIUpper) + ")");
//                buffer.append("{"); 
//                if (result.sasCILowerDeviation > sasTolerance)
//                    buffer.append("<font color='red'>" + Number.format(result.sasCILowerDeviation) + "</font>");
//                else
//                    buffer.append(Number.format(result.sasCILowerDeviation));
//                buffer.append(", ");
//                if (result.sasCIUpperDeviation > sasTolerance)
//                    buffer.append("<font color='red'>" + Number.format(result.sasCIUpperDeviation) + "</font>");
//                else
//                    buffer.append(Number.format(result.sasCIUpperDeviation));
//                buffer.append("}</td><td>");
//            }
//            buffer.append(Number.format(result.simulatedPower));
//            if (result.simulationDeviation > simTolerance)
//                buffer.append(" <font color='red'>(" + Number.format(result.simulationDeviation) + ")</font></td><td>");
//            else
//                buffer.append(" (" + Number.format(result.simulationDeviation) + ")</td><td>");
//
//            buffer.append(result.calculatedPower.getTest() + "</td><td>" + 
//                    Number.format(result.calculatedPower.getSigmaScale()) + "</td><td>" + 
//                    Number.format(result.calculatedPower.getBetaScale()) + "</td><td>" + 
//                    result.calculatedPower.getTotalSampleSize() + "</td><td>" + 
//                    Number.format(result.calculatedPower.getAlpha()) + "</td><td>" + 
//                    powerMethodToString(result.calculatedPower.getPowerMethod()) + "</td><td>" + 
//                    result.calculatedPower.getQuantile() + "</td></tr>");
//        }
            section.add(table);
        }

    }
    
}
