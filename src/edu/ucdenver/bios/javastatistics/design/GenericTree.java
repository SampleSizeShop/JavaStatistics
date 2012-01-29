/*
 * Study Design Service for the GLIMMPSE Software System.  
 * This service stores study design definitions for users of the GLIMMSE interface.
 * Service contain all information related to a power or sample size calculation.  
 * The Study Design Service simplifies communication between different screens in the user interface.
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
package edu.ucdenver.bios.javastatistics.design;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
 
public class GenericTree<T> {
 
    public T data;
    public List<GenericTree<T>> children;
 
    public GenericTree() {
        super();
        children = new ArrayList<GenericTree<T>>();
    }
 
    public GenericTree(T data) {
        this();
        setData(data);
    }
 
   public List<GenericTree<T>> getChildren() {
        return this.children;
    }
 
    public int getNumberOfChildren() {
        return getChildren().size();
    }
 
    public boolean hasChildren() {
        return (getNumberOfChildren() > 0);
    }
 
    public void setChildren(List<GenericTree<T>> children) {
        this.children = children;
    }
 
    public void addChild(GenericTree<T> child) {
        children.add(child);
    }
 
    public void addChildAt(int index, GenericTree<T> child) throws IndexOutOfBoundsException {
        children.add(index, child);
   }
 
    public void removeChildren() {
        this.children = new ArrayList<GenericTree<T>>();
   }
 
    public void removeChildAt(int index) throws IndexOutOfBoundsException {
        children.remove(index);
    }
 
    public GenericTree<T> getChildAt(int index) throws IndexOutOfBoundsException {
        return children.get(index);
    }
 
    public T getData() {
        return this.data;
    }
 
    public void setData(T data) {
        this.data = data;
    }
 
    public String toString() {
        return getData().toString();
    }
 
    public boolean equals(GenericTree<T> node) {
        return node.getData().equals(getData());
    }
 
    public int hashCode() {
        return getData().hashCode();
    }
 
    public String toStringVerbose() {
        String stringRepresentation = getData().toString() + ":[";
 
        for (GenericTree<T> node : getChildren()) {
            stringRepresentation += node.getData().toString() + ", ";
        }
 
        //Pattern.DOTALL causes ^ and $ to match. Otherwise it won't. It's retarded.
        Pattern pattern = Pattern.compile(", $", Pattern.DOTALL);
        Matcher matcher = pattern.matcher(stringRepresentation);
 
        stringRepresentation = matcher.replaceFirst("");
        stringRepresentation += "]";
 
        return stringRepresentation;
    }
}