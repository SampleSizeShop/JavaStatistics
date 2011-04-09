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
package edu.cudenver.bios.power.test.paper;

import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.xml.sax.InputSource;

import edu.cudenver.bios.distribution.ChiSquareTerm;
import edu.cudenver.bios.distribution.WeightedSumOfNoncentralChiSquaresDistribution;
import junit.framework.TestCase;

/**
 * Test case which compares the output of Davies Algorithm with 
 * the SAS/IML implementation of Davies.
 * 
 * @author Sarah Kreidler
 *
 */
public class TestDaviesAlgorithm extends TestCase
{
	private static final String TAG_TEST_CASE = "testCase";
	private static final String ATTR_ACCURACY = "accuracy";
	private static final String TAG_CHI_SQUARE = "chiSquare";
	private static final String ATTR_LAMBDA = "lambda";
	private static final String ATTR_OMEGA = "omega";
	private static final String ATTR_DF = "df";
	private static final String TAG_CDF_RESULT = "cdfResult";
	private static final String ATTR_QUANTILE = "quantile";
	private static final String ATTR_PROBABILITY = "probability";
	private static final String DATA_FILE = "data" + File.separator + "TestDaviesAlgorithm.xml";
    private static DecimalFormat Number = new DecimalFormat("#0.0000000");
    
	private ArrayList<ChiSquareTerm> chiSquareTerms = new ArrayList<ChiSquareTerm>();
	private HashMap<Double,Double> sasOutputValues = new HashMap<Double,Double>();
	private ArrayList<Double> quantileList = new ArrayList<Double>();
	private double accuracy = 0.001;
	
	WeightedSumOfNoncentralChiSquaresDistribution dist = null;
	
	/**
	 * Build the distribution and load the expected SAS values
	 */
	public void setUp()
	{
		FileReader reader = null;
        try
        {
            DocumentBuilderFactory factory =  DocumentBuilderFactory.newInstance();

            DocumentBuilder builder = factory.newDocumentBuilder();
            reader = new FileReader(DATA_FILE);
            Document doc = builder.parse(new InputSource(reader));

            NodeList testCaseNodes =  doc.getElementsByTagName(TAG_TEST_CASE);
            if (testCaseNodes.getLength() > 0)
            {
            	NamedNodeMap attrList = testCaseNodes.item(0).getAttributes();
            	Node accuracyNode = attrList.getNamedItem(ATTR_ACCURACY);
            	if (accuracyNode != null) accuracy = Double.parseDouble(accuracyNode.getNodeValue());
            }
            
            NodeList chiSquareNodes = doc.getElementsByTagName(TAG_CHI_SQUARE);
            for(int i = 0; i < chiSquareNodes.getLength(); i++)
            {
            	Node chiSquareNode = chiSquareNodes.item(i);
            	NamedNodeMap attrs = chiSquareNode.getAttributes();
            	Node lambdaNode = attrs.getNamedItem(ATTR_LAMBDA);
            	Node omegaNode = attrs.getNamedItem(ATTR_OMEGA);
            	Node dfNode = attrs.getNamedItem(ATTR_DF);
            	if (lambdaNode != null && omegaNode != null && dfNode != null)
            	{
            		chiSquareTerms.add(new ChiSquareTerm(Double.parseDouble(lambdaNode.getNodeValue()), 
            				Double.parseDouble(dfNode.getNodeValue()), Double.parseDouble(omegaNode.getNodeValue())));
            	}
            }
            
    		dist	= new WeightedSumOfNoncentralChiSquaresDistribution(chiSquareTerms, accuracy);

            NodeList resultNodes = doc.getElementsByTagName(TAG_CDF_RESULT);
            for(int i = 0; i < resultNodes.getLength(); i++)
            {
            	Node resultNode = resultNodes.item(i);
            	NamedNodeMap attrs = resultNode.getAttributes();
            	Node quantileNode = attrs.getNamedItem(ATTR_QUANTILE);
            	Node probabilityNode = attrs.getNamedItem(ATTR_PROBABILITY);
            	if (quantileNode != null && probabilityNode != null)
            	{
            		quantileList.add(Double.parseDouble(quantileNode.getNodeValue()));
            		sasOutputValues.put(Double.parseDouble(quantileNode.getNodeValue()), 
            				Double.parseDouble(probabilityNode.getNodeValue()));
            	}
            }
        }
        catch (Exception e)
        {
            fail();
        }
        finally
        {
        	if (reader != null) try { reader.close();} catch (Exception e) {};
        }
	}
	
	/**
	 * Calculate the probability for various quantiles for the given 
	 * distribution of chi-squares and compare against the SAS output
	 */
	public void testQuantile()
	{
		int matches = 0;
		for(double quantile: quantileList)
		{
			double p = dist.cdf(quantile);
			Double sasProb = sasOutputValues.get(quantile);
			System.out.println(quantile + "\t" + Number.format(p) + "\t(SAS value=" + sasProb + ")");
			if (!Double.isNaN(sasProb) && (Math.abs(sasProb - p) < accuracy))
				matches++;
		}
		
		assertEquals(quantileList.size(), matches);
	}
}
