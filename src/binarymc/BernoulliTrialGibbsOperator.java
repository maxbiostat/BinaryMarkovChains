package binarymc;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.util.Randomizer;


@Description("Bernoulli trial for randomly selected bit in an array.")
public class BernoulliTrialGibbsOperator extends Operator {
    final public Input<Double> pSuccessInput = new Input<>("p", "probability of success for a", 0.5);

    final public Input<BooleanParameter> parameterInput = new Input<>("parameter", "the parameter to operate a flip on.", Validate.REQUIRED);

    double pSuccess;
    BooleanParameter p;
    
    @Override
	public void initAndValidate() {
    	pSuccess = pSuccessInput.get();
    	p = parameterInput.get(this);
    }

    /**
     * Change the parameter and return the hastings ratio.
     * Assign value with probability p to a random bit in a bit vector.
     */

    @Override
    public double proposal() {

        final int dim = p.getDimension();

        final int pos = Randomizer.nextInt(dim);

        if (Randomizer.nextDouble() < pSuccess) {
            p.setValue(pos, true);
        } else {
            p.setValue(pos, false);
        }
        
        // make sure it is always accepted
        return Double.POSITIVE_INFINITY;
    }
}

