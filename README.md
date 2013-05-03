## Dynare for Julia

This package aims at bringing to Julia some of the algorithms incorporated in
[Dynare](http://www.dynare.org), a platform for solving dynamic economic
models.

Please note that this Julia package is very experimental and very incomplete
compared to the original Dynare for MATLAB/Octave. Use it only for experiment
purposes.

For the moment the package is only able to compute a model’s steady state and
first order decision rules. More features are expected in the future.

## A complete example

The following example is a translation of `example1.mod` that is shipped with
Dynare. It computes the steady state and first order decision rules.

Note that the syntax for the model block has been kept as close as possible to
the original format.

```
require("Dynare")

using Dynare

# Create a model object
m =
@modfile begin
    @var y c k a h b
    @varexo e u
    @parameters beta rho alpha delta theta psi tau
    @model begin
        c*theta*h^(1+psi) = (1-alpha)*y
        k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))*(exp(b(+1))*alpha*y(+1)+(1-delta)*k))
        y = exp(a)*(k(-1)^alpha)*(h^(1-alpha))
        k = exp(b)*(y-c)+(1-delta)*k(-1)
        a = rho*a(-1)+tau*b(-1) + e
        b = tau*a(-1)+rho*b(-1) + u
    end
end

# Do some preliminary computations
compute_model_info(m)

# Define a calibration and some starting values for the nonlinear solver

calib = [
         :alpha => 0.36,
         :rho   => 0.95,
         :tau   => 0.025,
         :beta  => 0.99,
         :delta => 0.025,
         :psi   => 0.0,
         :theta => 2.95
        ]

initval = [
           :y => 1.08068253095672,
           :c => 0.80359242014163,
           :h => 0.29175631001732,
           :k => 11.08360443260358,
           :a => 0.0,
           :b => 0.0,
           :e => 0.0,
           :u => 0.0
          ]

# Compute and print the steady state for the given calibration

s = steady_state(m, calib, initval)

print_steady_state(m, s)

# Compute and print eigenvalues and first order decision rules

(gy, gu, eigs) = decision_rules(m, calib, s)

println("Eigenvalues: ", eigs)
println()

print_decision_rules(m, gy, gu)
```
