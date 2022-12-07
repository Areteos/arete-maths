package calculus.differentialEquations;

import functions.MathsFunctions;
import types.tuples.Pair;

import java.util.*;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;


/**
 * Numerically solves partial differential equations of the form:<p>
 * &emsp&emsp  f<sub>t</sub> = Af<sub>xx</sub> + Bf<sub>x</sub> + Cf<p>
 * &emsp&emsp where f is an unknown function of t and x, and A, B, and C are known functions of t and x<p>
 * And returns a function that closely approximates the original function f.
 */
public abstract class PDESolver implements DoubleBinaryOperator {

    protected final BoundaryCondition lowerBound;
    protected final BoundaryCondition upperBound;
    protected final double spaceStep;
    protected final double timeStep;
    protected final double[] nodeLocations;
    protected final List<double[]> nodeSolutionsByTimeStep = new ArrayList<>();
    protected double currentTime = 0;
    protected final List<Double> timeSteps = new ArrayList<>();

    protected PDESolver(DoubleUnaryOperator functionAtT0, BoundaryCondition lowerBound, BoundaryCondition upperBound, double maxSpaceStep, double minTimeStep) {
        if (lowerBound.location() > upperBound.location()) {
            throw new IllegalArgumentException("Lower bound higher than upper bound");
        }

        this.lowerBound = lowerBound;
        this.upperBound = upperBound;

        double lowerBoundLocation = this.lowerBound.location();
        double upperBoundLocation = this.upperBound.location();

        double range = upperBoundLocation - lowerBoundLocation;
        int minimumSpaceSteps = (int) (range / maxSpaceStep) + 1;

        this.spaceStep = range / minimumSpaceSteps;
        nodeLocations = new double[minimumSpaceSteps+1];
        double[] initialNodeValues = new double[minimumSpaceSteps+1];
        nodeLocations[0] = lowerBoundLocation;
        initialNodeValues[0] = functionAtT0.applyAsDouble(nodeLocations[0]);
        for (int i = 1; i < minimumSpaceSteps; i++) {
            nodeLocations[i] = lowerBoundLocation;
            initialNodeValues[i] = functionAtT0.applyAsDouble(nodeLocations[i]);
        }

        nodeLocations[minimumSpaceSteps] = upperBoundLocation;
        initialNodeValues[minimumSpaceSteps] = functionAtT0.applyAsDouble(nodeLocations[minimumSpaceSteps]);

        nodeSolutionsByTimeStep.add(initialNodeValues);

        this.timeStep = minTimeStep;
        timeSteps.add(0d);
    }

    protected PDESolver(List<Pair<Double, Double>> samplesAndWeights, BoundaryCondition lowerBound, BoundaryCondition upperBound, double maxSpaceStep, double minTimeStep) {
        if (lowerBound.location() > upperBound.location()) {
            throw new IllegalArgumentException("Lower bound higher than upper bound");
        }

        List<Pair<Double, Double>> sortedSamplesAndWeights = new ArrayList<>(samplesAndWeights);
        sortedSamplesAndWeights.sort(Comparator.comparing(Pair::first));

        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        double lowerBoundLocation = this.lowerBound.location();
        double upperBoundLocation = this.upperBound.location();

        double range = upperBoundLocation - lowerBoundLocation;
        int minimumSpaceSteps = (int) (range / maxSpaceStep) + 1;
        this.spaceStep = range / minimumSpaceSteps;

        nodeLocations = new double[minimumSpaceSteps+1];
        double[] initialNodeValues = new double[minimumSpaceSteps+1];

        nodeLocations[0] = lowerBoundLocation;
        for (int i = 1; i < minimumSpaceSteps; i++) {
            nodeLocations[i] = nodeLocations[i-1] + spaceStep;
        }
        nodeLocations[minimumSpaceSteps] = upperBoundLocation;

        Iterator<Pair<Double, Double>> sampleAndWeightIterator = sortedSamplesAndWeights.iterator();
        Pair<Double, Double> sampleAndWeight = sampleAndWeightIterator.next();
        for (int i = 0; i <= minimumSpaceSteps && sampleAndWeight != null; i++) {
            double currentLocation = nodeLocations[i];

            while (MathsFunctions.isWithin(currentLocation, spaceStep, sampleAndWeight.first())) {
                initialNodeValues[i] += sampleAndWeight.second();
                if (sampleAndWeightIterator.hasNext()) {
                    sampleAndWeight = sampleAndWeightIterator.next();
                } else {
                    sampleAndWeight = null;
                    break;
                }
            }
        }

        nodeSolutionsByTimeStep.add(initialNodeValues);

        this.timeStep = minTimeStep;
        timeSteps.add(0d);
    }

    /**
     * Internally evaluates the solution for the PDE up to a given point in time. This function will be called internally
     * to lazily evaluate the solution on-demand.
     * @param time The time up to which to solve the equation
     */
    public abstract void solveUpToTime(double time);

    /**
     * Returns the value of the numerical solution to the PDE at a particular point in time and space.
     * @param t The time at which to find the solution
     * @param x The location at which to find the solution
     * @return The solution at the given time and location
     */
    @Override
    public double applyAsDouble(double t, double x) {
        if (t < 0) {
            throw new RuntimeException("Reverse timestep functionality not yet implemented");
        }

        if (t >= currentTime) {
            solveUpToTime(t);
        }
        int timeStepBeforeT = (int) (t/timeStep);
        int nodeBeforeX = (int) (x/spaceStep);

        if (nodeBeforeX >= nodeLocations.length-1) {
            nodeBeforeX = nodeLocations.length-2;
        }
        if (timeStepBeforeT >= timeSteps.size()-1) {
            timeStepBeforeT = timeSteps.size()-2;
        }
        if (nodeBeforeX < 0) {
            nodeBeforeX = 0;
        }
        if (timeStepBeforeT < 0) {
            timeStepBeforeT = 0;
        }

        double[] solutionsBeforeT = nodeSolutionsByTimeStep.get(timeStepBeforeT);
        double[] solutionsAfterT = nodeSolutionsByTimeStep.get(timeStepBeforeT + 1);
        double interpolatedSolutionBeforeT = MathsFunctions.interpolateLinearly(nodeLocations[nodeBeforeX], nodeLocations[nodeBeforeX+1], solutionsBeforeT[nodeBeforeX], solutionsBeforeT[nodeBeforeX+1], x);
        double interpolatedSolutionAfterT = MathsFunctions.interpolateLinearly(nodeLocations[nodeBeforeX], nodeLocations[nodeBeforeX+1], solutionsAfterT[nodeBeforeX], solutionsAfterT[nodeBeforeX+1], x);

        return MathsFunctions.interpolateLinearly(timeSteps.get(timeStepBeforeT), timeSteps.get(timeStepBeforeT+1), interpolatedSolutionBeforeT, interpolatedSolutionAfterT, t);
    }
}
