package binarymc;

import beast.util.Randomizer;

public class TSESS {
	public double pHat; // probability state = 1
	public double alphaHat; // probability of sticking at same state
	public double delta; // number of transitions
	public int sampleCount; // number of raw (correlated) samples

	public TSESS(Double [] trace, int burninPercentage) {
		int burnin = burninPercentage * trace.length / 100;
		sampleCount = trace.length - burnin;
		
		// estimate pHat
		pHat = 0;
		for (int i = burnin; i < trace.length; i++) {
			if (trace[i] != 0.0) {
				pHat++;
			}
		}
		pHat /= sampleCount;
		
		// estimate MC
		double p00=0,p01=0,p10=0,p11=0;
		for (int i = burnin; i < trace.length-1; i++) {
			if (trace[i] == 0.0) {
				if (trace[i+1] == 0.0) {
					p00++;
				} else {
					p01++;
				}
			} else {
				if (trace[i+1] == 0.0) {
					p10++;						
				} else {
					p11++;
				}
			}
		}
		// TODO: this is wrong -- plug in the right formula
		alphaHat = p01 * (1-pHat)/(p00+p01) + p10 * pHat/(p10+p11);
		delta = p01 + p10;
	}
	
	public double ESS() {
		// TODO: this is wrong -- plug in the right formula
		return sampleCount * alphaHat;
	}

	static Double [] resample(Double [] trace, int burninPercentage, int targetLength) {
		Double [] sample = new Double[targetLength];
		int [] length = new int[trace.length];
		for (int i = 0; i < targetLength; i++) {
			length[Randomizer.nextInt(trace.length)]++;
		}
		int k = 0;
		for (int i = 0; i < trace.length; i++) {
			Double v = trace[k];
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
} // class TSESS