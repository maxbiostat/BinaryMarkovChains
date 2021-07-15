package binarymc;

public class Test {

	private static double sqr(double x) {
		return x*x;
	}

	public static void main(String[] args) {
		double pHat = 0.5;
		
		// count Markov chain transitions
		double r = 1000; // repeat this experiment r times, so ESS = r * 4
		
		double  n00 = 99 * r, // stay 99 steps in state n00
				n01 = 1 * r,  // transition once
				n11 = n00,    // then stay 99 steps in state n11
				n10 = n01;    // transition once
	
		// k & v are shape parameters on the beta prior on alpha
		int k = 1;
		int v = 1;
		double z = (1-pHat)/pHat;
		double U = n00 + v - 1;
		double W = n01 + n10 + k - 1;
		// M = observation count
		double M = n00 + n01 + n10 + n11;

		double s = Math.sqrt(sqr((n11+W)*z) + 2*((U-W)*n11 - W*W - U*W)*z + sqr(U+W));
		
		double alpha = -(s + (-n11-W)*z-W-U)/(2 * M * z);
		
		double beta =alpha * (1-pHat)/pHat;

		System.out.println("a = " + alpha + " b = " + beta);
		System.out.println("M = " + M);
		System.out.println("estimated ESS = " + M*alpha/(2*pHat-alpha));
		System.out.println("correct   ESS = " + M*alpha/pHat);
//		double entropy = M*(-(1-alpha)*Math.log(1-alpha) - alpha*Math.log(alpha) +
//						 -(1-beta)*Math.log(1-beta) - beta*Math.log(beta));
//		double entropy2 = -p00*Math.log(p00/M)  -p01*Math.log(p01/M) 
//				 -p10*Math.log(p10/M)  -p11*Math.log(p11/M) ;
//		System.out.println("entropy = " + entropy + " " + entropy2);
	}
}
