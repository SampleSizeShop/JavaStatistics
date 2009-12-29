package edu.cudenver.bios.powersamplesize.graphics;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import edu.cudenver.bios.powersamplesize.Power;
import edu.cudenver.bios.powersamplesize.SampleSize;
import edu.cudenver.bios.powersamplesize.parameters.PowerSampleSizeParameters;

/**
 * 
 * @author kreidles
 *
 */
public class PowerCurveBuilder
{
	private static final double MAX_POWER = 0.999;
	private static final int MAX_DATA_POINTS = 10000; // safety net to avoid infinite loops on GLMM's
    // power and sample size calculator objects
    Power powerCalculator;
    SampleSize sampleSizeCalculator;
    // sample size increment and minimum - may be different than 1,2 for multivariate models
    int sampleSizeIncrement = 1;
    int minimumSampleSize = 2;
    // minimums and increment values for detectable difference (used if bySampleSize=false)
    double minimumMeanDifference = 0;
    double meanDifferenceIncrement = 1;
    // curve options (legend, axis labels, etc.)
    boolean legend = false;
    boolean bySampleSize = true;
    String title = "";
    String xaxisLabel = "";
    String yaxisLabel = "";
    String seriesName = "";
    // TODO
    
    /**
     * Constructor
     */
    public PowerCurveBuilder(Power powerCalculator, SampleSize sampleSizeCalculator) 
    {
        this.powerCalculator = powerCalculator;
        this.sampleSizeCalculator = sampleSizeCalculator;
    }

    public JFreeChart getPowerCurve(PowerSampleSizeParameters params)
    throws IllegalArgumentException
    {
    	// create a data series
    	XYSeries series = null;
    	if (bySampleSize)
    		series = buildSeriesPowerByN(params);
    	else
    		series = buildSeriesPowerbyMeanDifference(params);

        // complete the data set
        XYDataset powerData = new XYSeriesCollection(series);
        // use a spline renderer to make the connecting lines smooth
        XYSplineRenderer rend = new XYSplineRenderer();
        // turn off shapes displayed at each data point to make a smooth curve
        rend.setBaseShapesVisible(false);
        // Create the line chart
        XYPlot plot = new XYPlot(powerData, new NumberAxis(xaxisLabel), 
                new NumberAxis(yaxisLabel), rend);
        JFreeChart chart = new JFreeChart(title, 
                JFreeChart.DEFAULT_TITLE_FONT, plot, legend);

        return chart;
    }

    private XYSeries buildSeriesPowerByN(PowerSampleSizeParameters params)
    {
    	XYSeries series = new XYSeries(seriesName);
        /* determine end points of the power curve */
        // add a data point for power=alpha
        series.add(0, params.getAlpha());
        // get the sample size at which power nears 1 (>0.9999)
        params.setPower(MAX_POWER);
        int nMax = (int) Math.ceil(sampleSizeCalculator.getSampleSize(params));
        
        // calculate the power for various sample sizes and add them to the data series
        double calculatedPower = 0;
        for(int n = minimumSampleSize; n < nMax && calculatedPower < MAX_POWER; n += sampleSizeIncrement)
        {
            params.setSampleSize(n);
            calculatedPower = powerCalculator.getCalculatedPower(params);
            series.add(n, calculatedPower);
        }

        
        return series;
    }
    
    private XYSeries buildSeriesPowerbyMeanDifference(PowerSampleSizeParameters params)
    {
    	XYSeries series = new XYSeries(seriesName);
        /* determine end points of the power curve */
        // add a data point for power=alpha
        series.add(0, params.getAlpha());
        
        // calculate the power for various sample sizes and add them to the data series
        double calculatedPower = 0;
        int pts = 0;
        for(double delta = minimumMeanDifference; pts < MAX_DATA_POINTS && calculatedPower < MAX_POWER; pts++, delta += meanDifferenceIncrement)
        {
            params.setMeanDifference(delta);
            calculatedPower = powerCalculator.getCalculatedPower(params);
            series.add(delta, calculatedPower);
        }
        
        return series;    	
    }
    
	public int getSampleSizeIncrement()
	{
		return sampleSizeIncrement;
	}

	public void setSampleSizeIncrement(int sampleSizeIncrement)
	{
		this.sampleSizeIncrement = sampleSizeIncrement;
	}

	public int getMinimumSampleSize()
	{
		return minimumSampleSize;
	}

	public void setMinimumSampleSize(int minimumSampleSize)
	{
		this.minimumSampleSize = minimumSampleSize;
	}

	public boolean hasLegend()
	{
		return legend;
	}

	public double getMinimumMeanDifference()
	{
		return minimumMeanDifference;
	}

	public void setMinimumMeanDifference(double minimumMeanDifference)
	{
		this.minimumMeanDifference = minimumMeanDifference;
	}

	public double getMeanDifferenceIncrement()
	{
		return meanDifferenceIncrement;
	}

	public void setMeanDifferenceIncrement(double meanDifferenceIncrement)
	{
		this.meanDifferenceIncrement = meanDifferenceIncrement;
	}

	public void setLegend(boolean legend)
	{
		this.legend = legend;
	}

	public boolean isBySampleSize()
	{
		return bySampleSize;
	}

	public void setBySampleSize(boolean bySampleSize)
	{
		this.bySampleSize = bySampleSize;
	}

	public String getXaxisLabel()
	{
		return xaxisLabel;
	}

	public void setXaxisLabel(String xaxisLabel)
	{
		this.xaxisLabel = xaxisLabel;
	}

	public String getYaxisLabel()
	{
		return yaxisLabel;
	}

	public void setYaxisLabel(String yaxisLabel)
	{
		this.yaxisLabel = yaxisLabel;
	}

	public String getTitle()
	{
		return title;
	}

	public void setTitle(String title)
	{
		this.title = title;
	}
    
    
}
