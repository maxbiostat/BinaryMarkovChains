package binarymc;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.util.Randomizer;


@Description("Bernoulli trial for randomly selected bit in an array.")
public class BernoulliTrialGibbsOperator extends Operator {
    final public Input<Double> pSuccessInput = new Input<>("p", "probability of success for a", 0.5);
    final public Input<Boolean> randomOrderInput = new Input<>("randomOrder", "if true, randomly select a dimension, otherwise circle through dimensions deterministically", false);

    final public Input<BooleanParameter> parameterInput = new Input<>("parameter", "the parameter to operate a flip on.", Validate.REQUIRED);

    private double pSuccess;
    private boolean randomOrder;
    private boolean [] done;
    private int current = 0;
    private BooleanParameter p;
    
    @Override
	public void initAndValidate() {
    	pSuccess = pSuccessInput.get();
    	p = parameterInput.get(this);
    	randomOrder = randomOrderInput.get();
    	if (!randomOrder) {
    		done = new boolean[p.getDimension()];
    	}
    	
    }

    /**
     * Change the parameter and return the hastings ratio.
     * Assign value with probability p to a random bit in a bit vector.
     */

    @Override
    public double proposal() {

        final int dim = p.getDimension();

        final int pos = randomOrder ? 
        		Randomizer.nextInt(dim) :
        			next(dim);
//        			current++ % dim;

        if (Randomizer.nextDouble() < pSuccess) {
            p.setValue(pos, true);
        } else {
            p.setValue(pos, false);
        }
        
        // make sure it is always accepted
        return Double.POSITIVE_INFINITY;
    }

	private int next(int dim) {
		if (current == dim) {
			Arrays.fill(done, false);
			current = 0;
		}
		current++;
		int i = Randomizer.nextInt(dim);
		while (done[i]) {
			i++;
			if (i == dim) {
				i = 0;
			}
		}
		done[i] = true;
		return i;
	}
}

