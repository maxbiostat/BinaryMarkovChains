package binarymc;

import beast.util.Randomizer;

public class TSESS {
	public double pHat; // probability state = 1
	public double alphaHat; // probability of sticking at same state
	public double delta; // number of transitions
	public double M; // number of raw (correlated) samples

	// observation resample count for resampling as if resampleDeterministic() was called
	final static double r = 1;
	
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
		double n00=0,n01=0,n10=0,n11=0;
		for (int i = burnin; i < trace.length-1; i++) {
			if (trace[i] == 0.0) {
				if (trace[i+1] == 0.0) {
					n00 += 1;
				} else {
					n01 += 1;
				}
			} else {
				if (trace[i+1] == 0.0) {
					n10 += 1;						
				} else {
					n11 += 1;
				}
			}
		}
		
		// adjust constant observations as if the trace was 
		// resampled using the resampleDeterministic() method
		// with resampleCount = r
		n00 = n00 * r + n10 * (r-1);
		n11 = n11 * r + n01 * (r-1);
		
		// k & v are shape parameters on the beta prior on alpha
		int k = 1;
		int v = 1;
		double z = (1-pHat)/pHat;
		double U = n00 + v - 1;
		double W = n01 + n10 + k - 1;
		// M = observation count
		M = n00 + n01 + n10 + n11;

		double s = Math.sqrt(sqr((n11+W)*z) + 2*((U-W)*n11 - W*W - U*W)*z + sqr(U+W));
		
		alphaHat = -(s + (-n11-W)*z-W-U)/(2 * M * z);
		

		// alternative derivation of alphaHat, gives the same result
//		double a = (p00 + p01 + p10 + p11 + k - 1 + v - 1) * (1-pHat);
//		double b = - (U*pHat + W + p11*(1-pHat));
//		double c = W * pHat;
//		double alphaHat2 = (-b - Math.sqrt(b*b-4*a*c))/(2*a);
//		System.err.println(alphaHat + " " + (alphaHat - alphaHat2)/alphaHat);
		
		delta = n01 + n10;
	}
	
	private double sqr(double x) {
		return x*x;
	}
	
	public double ESS() {
//		if (r > 1) {
//			// multiply by fudge factor = 2
//			return 2 * M * alphaHat / (2*pHat -alphaHat);
//		}		
//		return M * alphaHat / (2*pHat -alphaHat);
		return M * alphaHat / pHat;
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