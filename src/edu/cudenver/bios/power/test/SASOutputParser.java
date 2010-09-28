package edu.cudenver.bios.power.test;

import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class SASOutputParser
{
    // statistical test name constants
    public static final String TEST_HOTELLING_LAWLEY_TRACE = "hlt";
    public static final String TEST_PILLAI_BARTLETT_TRACE = "pbt";
    public static final String TEST_WILKS_LAMBDA = "wl";
    public static final String TEST_UNIREP = "unirep";
    public static final String TEST_UNIREP_BOX = "unirepBox";
    public static final String TEST_UNIREP_GG = "unirepGG";
    public static final String TEST_UNIREP_HF = "unirepHF";
    
    // types of power approximations
    public static final String POWER_METHOD_CONDITIONAL = "conditional";
    public static final String POWER_METHOD_UNCONDITIONAL = "unconditional";
    public static final String POWER_METHOD_QUANTILE = "quantile";
    
	public static GLMMPowerParameters parseGLMMParameters(Document doc)
	{
		// TODO: have SAS code output the matrices it used to XML file
		return null;
	}
	
	public static List<GLMMPower> parsePowerResults(Document doc)
	{
		ArrayList<GLMMPower> powerList = new ArrayList<GLMMPower>();
		
		NodeList nodes = doc.getElementsByTagName("glmmPower");
		for(int i = 0; i < nodes.getLength(); i++)
		{
			Node glmmPower = nodes.item(i);
			NamedNodeMap attrs = glmmPower.getAttributes();

			Node testNode = attrs.getNamedItem("test");
			Node actualPowerNode = attrs.getNamedItem("actualPower");
			Node sampleSizeNode = attrs.getNamedItem("sampleSize");
			Node betaScaleNode = attrs.getNamedItem("betaScale");
			Node sigmaScaleNode = attrs.getNamedItem("sigmaScale");
			Node alphaNode = attrs.getNamedItem("alpha");
			Node nominalPowerNode = attrs.getNamedItem("nominalPower"); 
			Node powerMethodNode = attrs.getNamedItem("powerMethod");
			Node quantileNode = attrs.getNamedItem("quantile");

			if (testNode != null && actualPowerNode != null && 
					sampleSizeNode != null && betaScaleNode != null && 
					sigmaScaleNode != null && alphaNode != null &&
					nominalPowerNode != null && powerMethodNode != null)
			{
				if (quantileNode != null)
				{
					powerList.add(new GLMMPower(stringToTest(testNode.getNodeValue().trim()), 
							Double.parseDouble(alphaNode.getNodeValue().trim()), 
							Double.parseDouble(nominalPowerNode.getNodeValue().trim()), 
							Double.parseDouble(actualPowerNode.getNodeValue().trim()), 
							Integer.parseInt(sampleSizeNode.getNodeValue().trim()),
							Double.parseDouble(betaScaleNode.getNodeValue().trim()), 
							Double.parseDouble(sigmaScaleNode.getNodeValue().trim()), 
							stringToPowerMethod(powerMethodNode.getNodeValue().trim()),
							Double.parseDouble(quantileNode.getNodeValue().trim())));
				}
				else
				{
					powerList.add(new GLMMPower(stringToTest(testNode.getNodeValue().trim()), 
							Double.parseDouble(alphaNode.getNodeValue().trim()), 
							Double.parseDouble(nominalPowerNode.getNodeValue().trim()), 
							Double.parseDouble(actualPowerNode.getNodeValue().trim()), 
							Integer.parseInt(sampleSizeNode.getNodeValue().trim()),
							Double.parseDouble(betaScaleNode.getNodeValue().trim()), 
							Double.parseDouble(sigmaScaleNode.getNodeValue().trim()), 
							stringToPowerMethod(powerMethodNode.getNodeValue().trim())));
				}
			}
			else
			{
				System.err.println("Warning: invalid glmmPower object in input file");
			}
		}
		
		return powerList;
	}
	
	public static GLMMPowerParameters.PowerMethod stringToPowerMethod(String methodString)
	{
		if (methodString != null)
		{
			if (POWER_METHOD_QUANTILE.equals(methodString))
			{
    			return GLMMPowerParameters.PowerMethod.QUANTILE_POWER;
			}
			else if (POWER_METHOD_UNCONDITIONAL.equals(methodString))
			{
				return GLMMPowerParameters.PowerMethod.UNCONDITIONAL_POWER;
			}
			else if (POWER_METHOD_CONDITIONAL.equals(methodString))
			{
				return GLMMPowerParameters.PowerMethod.CONDITIONAL_POWER;
			}
		}
		return null;
	}
	
	public static GLMMPowerParameters.Test stringToTest(String testString)
	{
		if (testString != null)
		{
			if (TEST_HOTELLING_LAWLEY_TRACE.equals(testString))
			{
				return GLMMPowerParameters.Test.HOTELLING_LAWLEY_TRACE;
			}
			else if (TEST_PILLAI_BARTLETT_TRACE.equals(testString))
			{
				return GLMMPowerParameters.Test.PILLAI_BARTLETT_TRACE;
			}
			else if (TEST_WILKS_LAMBDA.equals(testString))
			{
				return GLMMPowerParameters.Test.WILKS_LAMBDA;
			}
			else if (TEST_UNIREP.equals(testString))
			{
				return GLMMPowerParameters.Test.UNIREP;
			}
			else if (TEST_UNIREP_BOX.equals(testString))
			{
				return GLMMPowerParameters.Test.UNIREP_BOX;
			}
			else if (TEST_UNIREP_GG.equals(testString))
			{
				return GLMMPowerParameters.Test.UNIREP_GEISSER_GREENHOUSE;
			}
			else if (TEST_UNIREP_HF.equals(testString))
			{
				return GLMMPowerParameters.Test.UNIREP_HUYNH_FELDT;
			}
		}
		return null;
	}
}
