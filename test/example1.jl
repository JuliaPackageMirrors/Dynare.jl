include("Dynare.jl")

m =
Dynare.@modfile begin
    Dynare.@var y c k a h b
    Dynare.@varexo e u
    Dynare.@parameters beta rho alpha delta theta psi tau
    Dynare.@model begin
        c*theta*h^(1+psi) = (1-alpha)*y
        k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))*(exp(b(+1))*alpha*y(+1)+(1-delta)*k))
        y = exp(a)*(k(-1)^alpha)*(h^(1-alpha))
        k = exp(b)*(y-c)+(1-delta)*k(-1)
        a = rho*a(-1)+tau*b(-1) + e
        b = tau*a(-1)+rho*b(-1) + u
    end
end

Dynare.compute_model_info(m)

println(m)

calib = [ :alpha => 0.36,
         :rho   => 0.95,
         :tau   => 0.025,
         :beta  => 0.99,
         :delta => 0.025,
         :psi   => 0.0,
         :theta => 2.95 ]

initval = [
           :y => 1.08068253095672,
           :c => 0.80359242014163,
           :h => 0.29175631001732,
           :k => 11.08360443260358,
           :a => 0.0,
           :b => 0.0,
           :e => 0.0,
           :u => 0.0 ]

s = Dynare.steady(m, calib, initval)

println("Steady state: ", s)

(gy, gu, eigs) = Dynare.decision_rules(m, calib, s)

println("Eigenvalues: ", eigs)

println(gy)

println(gu)
