package binarymc;

import beast.util.Randomizer;

public class TSESS {
	public double pHat; // probability state = 1
	public double alphaHat; // probability of sticking at same state
	public double delta; // number of transitions
	public double M; // number of raw (correlated) samples

	// observation resample count for resampling as if resampleDeterministic() was called
	final static double r = 1000;
	
	public TSESS(Double [] trace, int burninPercentage) {
		int burnin = burninPercentage * trace.length / 100;
		int sampleCount = trace.length - burnin;
		
		// estimate pHat
		pHat = 0;
		for (int i = burnin; i < trace.length; i++) {
			if (trace[i] != 0.0) {
				pHat++;
			}
		}
		pHat /= sampleCount;
		
		// count Markov chain transitions
		double p00=0,p01=0,p10=0,p11=0;
		for (int i = burnin; i < trace.length-1; i++) {
			if (trace[i] == 0.0) {
				if (trace[i+1] == 0.0) {
					p00 += 1;
				} else {
					p01 += 1;
				}
			} else {
				if (trace[i+1] == 0.0) {
					p10 += 1;						
				} else {
					p11 += 1;
				}
			}
		}
		
		// adjust constant observations as if the trace was 
		// resampled using the resampleDeterministic() method
		// with resampleCount = r
		p00 = p00 * r + p10 * (r-1);
		p11 = p11 * r + p01 * (r-1);
		
		// k & v are shape parameters on the beta prior on alpha
		int k = 1;
		int v = 1;
		double z = (1-pHat)/pHat;
		double U = p00 + v - 1;
		double W = p01 + p10 + k - 1;
		// M = observation count
		M = p00 + p01 + p10 + p11;

		double s = Math.sqrt(sqr((p11+W)*z) + 2*((U-W)*p11 - W*W - U*W)*z + sqr(U+W));
		
		alphaHat = -(s + (-p11-W)*z-W-U)/(2 * M * z);
		delta = p01 + p10;
	}
	
	private double sqr(double x) {
		return x*x;
	}
	
	public double ESS() {
		if (r > 1) {
			// multiply by fudge factor = 2
			return 2 * M * alphaHat / (2*pHat -alphaHat);
		}		
		return M * alphaHat / (2*pHat -alphaHat);
	}

	/**
	 * generate new trace with length equal to targetLength
	 * using the method from Luiza Fabreti & Sebastian Hoehna,
	 * https://www.biorxiv.org/content/10.1101/2021.05.04.442586v1
	 * @param trace input trace
	 * @param burninPercentage amount of burnin to disregards
	 * @param targetLength length of output trace
	 * @return trace with length targetLength resampled with replacement and sorted
	 */
	static Double [] resample(Double [] trace, int burninPercentage, int targetLength) {
		Double [] sample = new Double[targetLength];
		int [] length = new int[trace.length];
		for (int i = 0; i < targetLength; i++) {
			length[Randomizer.nextInt(trace.length)]++;
		}
		int k = 0;
		for (int i = 0; i < trace.length; i++) {
			Double v = trace[i];
			for (int j = 0; j < length[i]; j++) {
				sample[k] = v;
				k++;
			}
		}
		if (k != sample.length) {
			throw new RuntimeException("Programmer error: k should be equal to target length");
		}
		return sample;
	}
	
	/**
	 * generate new trace where each observation is multiplied resampleCount times
	 * can be considered a deterministic version of resample() above
	 */
	static Double [] resampleDeterministic(Double [] trace, int burninPercentage, int resampleCount) {
		Double [] sample = new Double[resampleCount * trace.length];
		int k = 0;
		for (int i = 0; i < trace.length; i++) {
			Double v = trace[i];
			for (int j = 0; j < resampleCount; j++) {
				sample[k] = v;
				k++;
			}
		}
		if (k != sample.length) {
			throw new RuntimeException("Programmer error: k should be equal to target length");
		}
		return sample;
	}
	
} // class TSESS