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
package edu.ucdenver.bios.power;

import edu.cudenver.bios.power.GLMMPower;
import edu.cudenver.bios.power.test.SASOutputParser;
import org.w3c.dom.Document;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

public class Utils {

    /**
     * Read the XML file identified by resource name from the classpath, parse
     * and return a list of {@link GLMMPower}s.  Test fails if any exceptions
     * are thrown reading the XML file.
     *
     * @param resourceName name of resource to parse
     * @return List of {@link GLMMPower}s from XML file
     * @throws ParserConfigurationException
     * @throws IOException
     * @throws SAXException
     */
    public static List<GLMMPower> readSasPowers(String resourceName) {
        try {
            assertNotNull("resourceName cannot be null", resourceName);

            // parse the sas xml file
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            InputStream is = Utils.class.getClassLoader().getResourceAsStream(resourceName);
            assertNotNull(resourceName + " not found on classpath", is);
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(new InputSource(is));
            return SASOutputParser.parsePowerResults(doc);
        } catch (Exception e) {
            e.printStackTrace();
            fail(e.getMessage());
            return null;  // for the compiler's sake
        }
    }
}
