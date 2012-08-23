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
import java.io.IOException;

import com.itextpdf.text.BaseColor;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Font;
import com.itextpdf.text.pdf.PdfWriter;

/**
 * Utility routines for creating PDF documents via iText.
 * @author Sarah Kreidler
 *
 */
public class PDFUtils {
    public static Font SECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 18,
            Font.BOLD);
    public static Font SUBSECTION_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 16,
            Font.BOLD);
    public static Font RED_FONT= new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.NORMAL, BaseColor.RED);
    public static Font BOLD_FONT = new Font(Font.FontFamily.TIMES_ROMAN, 12,
            Font.BOLD);
    
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
    
    
}
