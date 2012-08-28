package edu.cudenver.bios.power.test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;

import org.xhtmlrenderer.pdf.ITextRenderer;

import com.lowagie.text.DocumentException;

public class HtmlToPdf {

    /**
     * Convert all html files in a directory to pdf.
     * 
     * @param args command line args
     */
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("usage: HtmlToPdf <directory-of-html-file>");
        } else {
            File dir = new File(args[0]);
            String[] children = dir.list();
            if (children == null) {
                // Either dir does not exist or is not a directory
                System.err.println("Directory does not exist or is empty");
            } else {
                for (int i=0; i<children.length; i++) {
                    // Get filename of file or directory
                    try {
                        // make sure we are converting an HTML file
                        String filename = dir.getAbsolutePath() + 
                                File.separatorChar + children[i];
                        if (filename.endsWith(".html")) {
                            convertToPdf(filename, filename.replace(".html", ".pdf"));
                        }
                    } catch (Exception e) {
                        System.err.println(e.getMessage());
                    }
                }
            }
        }
    }
    
    private static void convertToPdf(String inFile, String outFile) 
    throws DocumentException, MalformedURLException,
    FileNotFoundException {
        String url = new File(inFile).toURI().toURL().toString();
        OutputStream os = new FileOutputStream(outFile);
        
        ITextRenderer renderer = new ITextRenderer();
        renderer.setDocument(url);
        renderer.layout();
        renderer.createPDF(os);
    }

//    /**
//     * Convert the specified HTML file to a pdf
//     * @param filename
//     * @throws DocumentException
//     */
//    private static void convertToPdf(String filename) 
//    throws DocumentException {
//        // TODO Auto-generated method stub
//        Document document = new Document();
//        try {
//            System.out.println("Converting file ["+filename+"]");
//
//            // create the pdf file of the same name
//            String pdfFilename = filename.replace(".html", ".pdf");
//            PdfWriter writer = PdfWriter.getInstance(document,
//                    new FileOutputStream(pdfFilename));
//            writer.setInitialLeading(12.5f);
//            document.open();
//
//            // CSS
//            CSSResolver cssResolver = new StyleAttrCSSResolver();
//            // HTML
//            XMLWorkerFontProvider fontProvider = new XMLWorkerFontProvider();
//            CssAppliers cssAppliers = new CssAppliersImpl(fontProvider);
//            HtmlPipelineContext htmlContext = new HtmlPipelineContext(cssAppliers);
//            htmlContext.setTagFactory(Tags.getHtmlTagProcessorFactory());
//           
//            // Pipelines
//            PdfWriterPipeline pdf = new PdfWriterPipeline(document, writer);
//            HtmlPipeline html = new HtmlPipeline(htmlContext, pdf);
//            CssResolverPipeline css = new CssResolverPipeline(cssResolver, html);
//           
//            XMLWorker worker = new XMLWorker(css, true);
//            XMLParser p = new XMLParser(worker);
//            p.parse(new FileInputStream(filename));
//           
//            // step 5
//            document.close();
//
//            System.out.println("Done.");
//        } catch (Exception e) {
//            document.close();
//            throw new DocumentException("failed to convert file [" +
//                    filename + "]: " + e.getMessage());
//        }
//    }
}
