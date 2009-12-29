package edu.cudenver.bios.matrix;

/**
 * Container class for meta data associated with a matrix row.
 * Primarily used with a design "essence" matrix to specify repeated
 * rows for producing a full design matrix.
 * <p>
 * Repetitions can be specified in two ways:
 * <ul>
 * <li>Actual value - number of times this row is repeated in the full
 * design matrix.  Required for power calculations</li>
 * <li>Ratio value - ratio of group sizes in the full design matrix.  Required
 * for sample size calculations</li>
 * </ul>
 * 
 * @see EssenceMatrix
 * @author kreidles
 *
 */
public class RowMetaData
{
    int repetitions = -1;
    
    int ratio = 1;
    /**
     * Constructor, creates an empty row meta data object
     */
    public RowMetaData() {}
    
    /**
     * Constructor.  Create a RowMetaData object with the specified
     * number of row repetitions.
     * 
     * @param repetitions either actual value or ratio of repetitions of row in 
     * full design matrix.
     * @throws IllegalArgumentException
     */
    public RowMetaData (int repetitions, int ratio) 
    throws IllegalArgumentException
    {
        if (repetitions <= 0 || ratio <= 0) 
            throw new IllegalArgumentException("Repetitions must be greater than zero");
        this.repetitions = repetitions;
        this.ratio = ratio;
    }

    /**
     * Get the number of times this row is repeated in the full matrix
     * 
     * @return repetitions
     */
    public int getRepetitions()
    {
        return repetitions;
    }

    /**
     * Set the number of times this row is repeated repetitions for the row
     * 
     * @param repetitions
     */
    public void setRepetitions(int repetitions)
    throws IllegalArgumentException
    {
        if (repetitions <= 0) 
            throw new IllegalArgumentException("Repetitions must be greater than zero");
        this.repetitions = repetitions;
    }

	public int getRatio()
	{
		return ratio;
	}

	public void setRatio(int ratio)
	{
		this.ratio = ratio;
	}   
    
    
}
