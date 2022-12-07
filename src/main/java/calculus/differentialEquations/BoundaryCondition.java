/*
 *    Copyright 2022 Glenn Mamacos
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package calculus.differentialEquations;

import java.util.function.DoubleUnaryOperator;

/**
 * A record to encapsulate all salient information about a 1D boundary condition for a partial differential equation
 * @param boundaryType The type of condition being imposed on the boundary
 * @param location The spatial location of the boundary
 * @param value The value of the boundary condition, as a function of time
 * @param order0Coefficient The coefficient of the solution term of the boundary condition, as a function of time
 * @param order1Coefficient The coefficient of the derivative term of the boundary condition, as a function of time
 */
public record BoundaryCondition(BoundaryType boundaryType, double location, DoubleUnaryOperator value, DoubleUnaryOperator order0Coefficient, DoubleUnaryOperator order1Coefficient) {
    public enum BoundaryType {
        /**
         * Specifies a value that the solution must take at a boundary.
         */
        DIRICHLET,
        /**
         * Specifies a value that the derivative must take at a boundary.
         */
        NEUMANN,
        /**
         * Specifies a value that a particular linear combination of the solution and the derivative must take at a boundary.
         */
        ROBIN
    }

    public double value(double time) {
        return value.applyAsDouble(time);
    }

    public double order0Coefficient(double time) {
        return order0Coefficient.applyAsDouble(time);
    }

    public double order1Coefficient(double time) {
        return order1Coefficient.applyAsDouble(time);
    }

    /**
     * Constructs a Dirichlet boundary condition.
     * @param location The spatial location of the boundary
     * @param value The value of the boundary condition, as a function of time
     * @return A Dirichlet boundary condition
     * @see BoundaryType#DIRICHLET
     */
    public static BoundaryCondition dirichlet(double location, DoubleUnaryOperator value) {
        return new BoundaryCondition(BoundaryType.DIRICHLET, location, value, null, null);
    }

    /**
     * Constructs a Neumann boundary condition.
     * @param location The spatial location of the boundary
     * @param value The value of the boundary condition, as a function of time
     * @return A Dirichlet boundary condition
     * @see BoundaryType#NEUMANN
     */
    public static BoundaryCondition neumann(double location, DoubleUnaryOperator value) {
        return new BoundaryCondition(BoundaryType.NEUMANN, location, value, null, null);
    }

    /**
     * Constructs a Robin boundary condition.
     * @param location The spatial location of the boundary
     * @param value The value of the boundary condition, as a function of time
     * @param order0Coefficient The coefficient of the solution term of the boundary condition, as a function of time
     * @param order1Coefficient The coefficient of the derivative term of the boundary condition, as a function of time
     * @see BoundaryType#ROBIN
     */
    public static BoundaryCondition robin(double location, DoubleUnaryOperator value, DoubleUnaryOperator order0Coefficient, DoubleUnaryOperator order1Coefficient) {
        return new BoundaryCondition(BoundaryType.ROBIN, location, value, order0Coefficient, order1Coefficient);
    }
}
