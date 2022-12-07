package kernelDensityEstimation;

import functions.IterableFunctions;
import functions.MathsFunctions;
import calculus.differentialEquations.BoundaryCondition;
import calculus.differentialEquations.PDESolverVTCS;
import calculus.differentiation.DifferentiableFunction;
import jdk.jfr.Experimental;
import types.tuples.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

import static java.lang.Math.*;
import static functions.IterableFunctions.zipToList;

@Experimental
public class KDEViaDiffusion {
    private final double
            spaceFineness = 0.01,
            timeFineness = 0.001;

    private final double cumulativeWeight;
    final static double xi = pow((6 * sqrt(2) - 3) / 7, 0.4);
    final static double machinePrecision = 1e-10;
    private final List<Pair<Double, Double>> samplesAndWeights;
    private final double maxSpaceStep;
    private final double spaceStep;

    private KDEViaDiffusion(List<Pair<Double, Double>> samplesAndWeights) {
        this.samplesAndWeights = List.copyOf(samplesAndWeights);
        this.maxSpaceStep = spaceFineness / this.samplesAndWeights.size();
        int minimumSpaceSteps = (int) (1 / maxSpaceStep) + 1;
        this.spaceStep = 1d / minimumSpaceSteps;

        Pair<Double, Double> minMax = IterableFunctions.getMinimumAndMaximum(this.samplesAndWeights.stream().map(Pair::first).toList()).orElseThrow();
        if (minMax.first() < 0 || minMax.second() > 1) {
            throw new RuntimeException("Samples for diffusion estimator must be unit-normalised");
        }
        cumulativeWeight = IterableFunctions.getSum(this.samplesAndWeights.stream().map(Pair::second).toList());
    }

    /**
     * @param l Number of stages to use in the bandwidth estimation process. 5 Is usually optimal.
     * @return A pair: firstly the gaussian probability density function with optimum bandwidth,
     * secondly the square of the optimal intermediate bandwidth for selecting bandwidth for the diffusion estimator
     */
    public Pair<DifferentiableFunction, Double> improvedSleatherJones_Algorithm1(int l) {
        List<Double> z = new ArrayList<>();
        z.add(0, machinePrecision);
        double change = Double.MAX_VALUE;
        while (change >= machinePrecision || z.size() < 10) {
            z.add(xi * gamma(1, l, z.get(z.size()-1)));
            change = Math.abs(z.get(z.size()-1) - z.get(z.size()-2));
            System.out.printf("Delta: %s%n", change);
        }
        double squaredBandwidth = z.get(z.size()-1);

        DifferentiableFunction kernel = Kernels.gaussian(sqrt(squaredBandwidth));
        DifferentiableFunction kernelDensityEstimator = KernelDensityEstimator.getKernelDensityEstimator(samplesAndWeights, kernel);
        double squaredIntermediateBandwidth = gamma(1, l-1, squaredBandwidth);

        System.out.printf("Gaussian squared bandwidth: %s%n", squaredBandwidth);
        System.out.printf("Gaussian bandwidth: %s%n", sqrt(squaredBandwidth));

        return new Pair<>(kernelDensityEstimator, squaredIntermediateBandwidth);
    }


    public double squareOfNormOfJthDerivativeOfFForBandwidthJ_26(double squaredBandwidth, int j) {
        DifferentiableFunction gaussian = Kernels.gaussian(sqrt(2 * squaredBandwidth));
        DifferentiableFunction gaussianDerivative = gaussian.differentiate(j * 2);

        DifferentiableFunction factorisedGaussianDerivative = gaussianDerivative.factorise();

        double sum = 0;
        for (Pair<Double, Double> zipX : samplesAndWeights) {
            double sampleX = zipX.first();
            double weightX = zipX.second();
            for (Pair<Double, Double> zipY : samplesAndWeights) {
                double sampleY = zipY.first();
                double weightY = zipY.second();
                sum += weightX * weightY * factorisedGaussianDerivative.applyAsDouble(sampleX - sampleY);
            }
        }
        return sum * pow(-1, j) / pow(cumulativeWeight, 2);
    }

    public double estimateSquaredBandwidthJ_29(double squareBandwidthJPlus1, int j) {
        double leftFraction = (1 + 1 / pow(2, j+0.5))/3;

        double rightNumerator = 1;
        for (int i = 3; i <= 2*j - 1; i+=2) {
            rightNumerator *= i;
        }
        double rightDenominator = cumulativeWeight * sqrt(PI/2) * squareOfNormOfJthDerivativeOfFForBandwidthJ_26(squareBandwidthJPlus1, j+1);
        double rightFraction = rightNumerator / rightDenominator;

        double exponent = 2d / (3 + 2*j);

        return pow(leftFraction * rightFraction, exponent);
    }


    public double gamma(int j, int l, double squaredBandwidthLPLus1) {
        double currentSquaredBandwidth = squaredBandwidthLPLus1;
        for (int i = l; i >= j; i--) {
            currentSquaredBandwidth = estimateSquaredBandwidthJ_29(currentSquaredBandwidth, i);
        }
        return currentSquaredBandwidth;
    }



    // =================================================================================================================
    // DIFFUSION
    // =================================================================================================================


    /**
     * INCOMPLETE, DO NOT USE
     */
    public DoubleUnaryOperator getDiffusionEstimator_Algorithm2(double alpha, int l) {
        Pair<DifferentiableFunction, Double> gaussianPDFAndSquaredIntermediateBandwidth = improvedSleatherJones_Algorithm1(l);
        DifferentiableFunction gaussianPDF = gaussianPDFAndSquaredIntermediateBandwidth.first();
        double squaredIntermediateBandwidth = gaussianPDFAndSquaredIntermediateBandwidth.second();

        System.out.printf("Intermediate squared bandwidth: %s%n", squaredIntermediateBandwidth);
        System.out.printf("Intermediate bandwidth: %s%n", sqrt(squaredIntermediateBandwidth));

        double minTimeStep = timeFineness * squaredIntermediateBandwidth;


        DoubleUnaryOperator distributionFunction = x -> {
            double pointWeight = 0;
            for (Pair<Double, Double> sampleAndWeight : samplesAndWeights) {
                double sample = sampleAndWeight.first();
                double weight = sampleAndWeight.second();

                double diff = x - sample;

                if ((diff >= 0 && diff < spaceStep/2) || (diff <= 0 && diff >= -spaceStep/2)) {
                    pointWeight += weight;
                }
            }
            return pointWeight / (cumulativeWeight * spaceStep);
        };

        Map<Double, Double> aSolutions = new HashMap<>();
        Map<Double, Double> bSolutions = new HashMap<>();
        Map<Double, Double> cSolutions = new HashMap<>();

        DoubleBinaryOperator A = (t, x) -> aSolutions.computeIfAbsent(x, x0-> pow(gaussianPDF.applyAsDouble(x0), alpha - 1) / 2);

        DoubleBinaryOperator B = (t, x) -> bSolutions.computeIfAbsent(x, x0 -> gaussianPDF.applyAsDouble(x0) == 0 ? 0 :
                (alpha - 2) * gaussianPDF.differentiate().applyAsDouble(x0) * pow(gaussianPDF.applyAsDouble(x0), alpha - 2) / 2);

        DoubleBinaryOperator C = (t, x) -> cSolutions.computeIfAbsent(x, x0 -> gaussianPDF.applyAsDouble(x0) == 0 ? 0 :
                pow(gaussianPDF.applyAsDouble(x0), alpha - 2)
                * (-gaussianPDF.differentiate(2).applyAsDouble(x0)
                - (alpha - 2) * pow(gaussianPDF.differentiate().applyAsDouble(x0), 2) / gaussianPDF.applyAsDouble(x0))
                / 2);

        double gaussianDerivativeAt0 = gaussianPDF.differentiate().applyAsDouble(0);
        double gaussianDerivativeAt1 = gaussianPDF.differentiate().applyAsDouble(1);
        double gaussianAt0 = gaussianPDF.applyAsDouble(0);
        double gaussianAt1 = gaussianPDF.applyAsDouble(1);

        double a0 = gaussianDerivativeAt0 / pow(gaussianAt0, 2);
        double b0 = 1 / gaussianAt0;

        double a1 = gaussianDerivativeAt1 / pow(gaussianAt1, 2);
        double b1 = 1 / gaussianAt1;

        BoundaryCondition lowerBound = BoundaryCondition.robin(0, t->0, t->a0, t->-b0);
        BoundaryCondition upperBound = BoundaryCondition.robin(1, t->0, t->a1, t->-b1);

        DoubleBinaryOperator numericalSolutionForDiffusionEstimator = new PDESolverVTCS(1,
            distributionFunction, A, B, C, lowerBound, upperBound,
            maxSpaceStep, minTimeStep
        );

        double pointWiseBiasTerm = 0;
        for (Pair<Double, Double> sampleAndWeight : samplesAndWeights) {
            double sample = sampleAndWeight.first();
            double weight = sampleAndWeight.second();

            pointWiseBiasTerm += weight * pow(gaussianPDF.applyAsDouble(sample), (1 - alpha) / 2);
        }
        pointWiseBiasTerm /= cumulativeWeight;


        double squaredBandwidth = getSquareOfAsymptoticallyOptimalDiffusionBandwidth_23(
                approximateNormLfSquared(numericalSolutionForDiffusionEstimator, 0.001, squaredIntermediateBandwidth),
                pointWiseBiasTerm
        );

        if (!Double.isFinite(squaredBandwidth)) {
            throw new RuntimeException("Bandwidth cannot be NaN");
        }

        System.out.printf("Chosen squared bandwidth: %s%n", squaredBandwidth);
        System.out.printf("Chosen bandwidth: %s%n", sqrt(squaredBandwidth));

        return x -> numericalSolutionForDiffusionEstimator.applyAsDouble(squaredBandwidth, x);
    }

    public double getSquareOfAsymptoticallyOptimalDiffusionBandwidth_23(double normLfSquared, double pointWiseBiasTerm) {
        return pow(pointWiseBiasTerm / (2 * cumulativeWeight * Math.sqrt(PI) * normLfSquared), 2/5d);
    }

//    public DoubleUnaryOperator getSquareOfAsymptoticallyOptimalDiffusionBandwidth_23(DoubleUnaryOperator normLfSquared, DoubleUnaryOperator pointWiseBiasTerm) {
//        return x -> pow(pointWiseBiasTerm.applyAsDouble(x) / (2 * cumulativeWeight * Math.sqrt(PI) * normLfSquared.applyAsDouble(x)), 2/5d);
//    }


    public double approximateNormLfSquared(DoubleBinaryOperator numericalSolutionForDiffusionEstimator, double epsilon, double t) {
        if (epsilon >= 1) {
            throw new RuntimeException("Epsilon value should be significantly less than 1");
        }

        double sum = 0;
        for (Pair<Double, Double> sampleAndWeight : samplesAndWeights) {
            double sample = sampleAndWeight.first();
            double weight = sampleAndWeight.second();

            sum += weight * Math.pow(
                    numericalSolutionForDiffusionEstimator.applyAsDouble(t + epsilon, sample)
                    - numericalSolutionForDiffusionEstimator.applyAsDouble(t, sample), 2);
        }
        return pow(sum / epsilon, 2);
    }

    /**
     * INCOMPLETE, DO NOT USE
     */
    public static DoubleUnaryOperator getDiffusionEstimator(List<Double> samples, List<Double> weights) {
        Pair<Double, Double> minMax = IterableFunctions.getMinimumAndMaximum(samples).orElseThrow();
        return getDiffusionEstimator(samples, weights, minMax.first(), minMax.second());
    }

    /**
     * INCOMPLETE, DO NOT USE
     */
    public static DoubleUnaryOperator getDiffusionEstimator(List<Double> samples, List<Double> weights, double minimum, double maximum) {
        DoubleUnaryOperator normalisationFunction = MathsFunctions.getLinearMappingFunction(minimum, maximum, 0, 1).orElseThrow();
        double range = maximum - minimum;

        KDEViaDiffusion kdeViaDiffusion = new KDEViaDiffusion(zipToList(samples.stream().map(normalisationFunction::applyAsDouble).toList(), weights));

        DoubleUnaryOperator diffusionEstimator = kdeViaDiffusion.getDiffusionEstimator_Algorithm2(1, 5);

        return diffusionEstimator.compose(normalisationFunction).andThen(x -> x / range);
    }

    /**
     * Given a set of samples and their respective weights, create a Gaussian Kernel Density Estimator with an
     * automatically selected, asymptotically optimal bandwidth. This is done using the Improved Sleather Jones method.
     * @param samples The sample data from the probability distribution to be approximated
     * @param weights The relative weights of those data
     * @param lowerBound The lower bound of the probability distribution
     * @param upperBound The upper bound of the probability distribution
     * @return A kernel density estimator that approximates the probability density function that created the input samples
     */
    public static DoubleUnaryOperator getGaussianEstimator(List<Double> samples, List<Double> weights, double lowerBound, double upperBound) {
        DoubleUnaryOperator normalisationFunction = MathsFunctions.getLinearMappingFunction(lowerBound, upperBound, 0, 1).orElseThrow();
        double range = upperBound - lowerBound;

        KDEViaDiffusion kdeViaDiffusion = new KDEViaDiffusion(zipToList(samples.stream().map(normalisationFunction::applyAsDouble).toList(), weights));

        DifferentiableFunction estimator = kdeViaDiffusion.improvedSleatherJones_Algorithm1(5).first();

        return estimator.compose(normalisationFunction).andThen(x -> x / range);
    }
}
