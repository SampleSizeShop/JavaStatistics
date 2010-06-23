package edu.cudenver.bios.power;

public class PowerException extends Exception
{
	static final long serialVersionUID = -1L;
	
	public PowerException(String msg)
	{
		super(msg);
	}
	
	public PowerException(Throwable e)
	{
		super(e);
	}
}
