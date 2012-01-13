proc iml;

ndf = 4;
ddf = 5;
p = {0.05 0.25 0.5 0.75 0.95};

foo = FNONCT(1,2,5.6,0.95);
print foo;


do f = 1 to 4 by 1;
	do nc = 1 to 5 by 0.5;
		prob = probf(f,ndf,ddf,nc);
		nc = FNONCT(f, ndf, ddf, );
			HOLDNC = HOLDNC//(f || ndf[i] || ddf[j] || p[k] || nc);
		end;
	end;
end;

names = {"F", "NDF", "DDF", "P", "NC"};
print HOLDNC[COLNAME=names];


    	System.out.println("Testing noncentrality:");
    	int failures = 0;
    	double[] fList = {1,2,3,4};
    	double ndf = 4;
    	double ddf = 5;

    	for(double f : fList)
    	{
    		for(double nc = 1; nc <= 5; nc += 0.5)
    		{
    			NonCentralFDistribution dist = new NonCentralFDistribution(ndf, ddf, nc);
    			double prob = dist.cdf(f);
    			try
    			{
    				double newNC = NonCentralFDistribution.noncentrality(f, prob, ndf, ddf, 1);
    				System.out.println("Ndf=" + ndf + " ddf=" + ddf + " f=" + f + " prob=" + prob + " true nc=" + nc + " calc nc = " + newNC);
    				if (Math.abs(nc - newNC) > 1.0E-6) failures++;
    			}
    			catch (Exception e)
    			{
    				System.err.println("Noncentrality Failed: " + e.getMessage());
    				failures++;
    			}
    		}
    	}
    	assertEquals(0, failures);	



quit;
