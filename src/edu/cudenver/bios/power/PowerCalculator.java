package edu.cudenver.bios.power;

import java.util.List;

import edu.cudenver.bios.power.parameters.PowerParameters;

public interface PowerCalculator
{
	List<Power> getPower(PowerParameters params);
	
	List<Power> getSampleSize(PowerParameters params);
	
	List<Power> getDetectableDifference(PowerParameters params);
	
	List<Power> getSimulatedPower(PowerParameters params, int iterations);
}
