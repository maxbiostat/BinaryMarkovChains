package binarymc;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.List;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.util.LogAnalyser;

@Description("Produces two-state effective sample size estimates on binary trace logs")
public class TwoStateESSApp extends beast.core.Runnable {
	final public Input<File> traceInput = new Input<>("in", "trace file containing binary entries", new File("[[none]]"));
	public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified", new OutFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trace to be used as burn-in (and will be ignored)", 10);
	final public Input<Integer> resampleInput = new Input<>("resample", "number of times the chain should be resampled", 1);

	
	@Override
	public void initAndValidate() {
	}

	final static String space = "                                 ";
	
	@Override
	public void run() throws Exception {
		LogAnalyser traces = new LogAnalyser(traceInput.get().getPath(), burnInPercentageInput.get(), true, false);
		
		int resample = resampleInput.get();
		
        PrintStream out = System.out;
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
        	out = new PrintStream(outputInput.get());
        }
        
        List<String> labels = traces.getLabels();
        DecimalFormat f = new DecimalFormat("#.##");
    	out.println("label" + space.substring(5) +
    			"pHat" + "\t" +
    			"TSESS");
        for (int i = 0; i < labels.size(); i++) {
        	String label = labels.get(i);
        	Double [] trace = traces.getTrace(label);
        	if (resample > 1) {
//        		trace = TSESS.resampleDeterministic(trace, 0, resample);
        		trace = TSESS.resample(trace, 0, trace.length * resample);
        	}
        	TSESS tsESS = new TSESS(trace, 0);
        	out.print(label + space.substring(label.length()) +
        			f.format(tsESS.pHat) + "\t" +
        			f.format(tsESS.ESS()));
        	out.println();
        }
        
        
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
        	out.close();
        }
        Log.warning.println("Done.");
	}


	public static void main(String[] args) throws Exception {
		new Application(new TwoStateESSApp(), "Two State ESS", args);
	}
}
