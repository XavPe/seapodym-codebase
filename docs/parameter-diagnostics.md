# Parameter diagnostics in Hessian output

`Estimate` is a physical parameter value. The Hessian, its inverse, `Gradient`,
`StdErr`, `NewtonStep` and `NewtonStep/SE` use the optimizer's internal coordinates.

For a separate transformation theta_i = f_i(x_i), the physical standard error
is approximated by `abs(d theta_i / d x_i) * StdErr_internal_i`. The reported
CV is this physical standard error divided by `abs(Estimate_i)`.

The Jacobian uses ADMB's bound transformation derivative and the effective
parameter bounds. Estimates and Jacobians are captured at the central vector,
before Hessian perturbations. The output keeps its seven-column layout and
adds the metadata line `# CV_COORDINATES: PHYSICAL_DELTA_METHOD`.

A zero estimate or invalid standard error gives an undefined CV (`NaN`).
Inverse-Hessian uncertainty still requires a valid local covariance and an
appropriate stationary point. Linear propagation can degenerate at bounds.
Changing diagnostic coordinates does not validate an indefinite Hessian.

The computation is confined to diagnostics and does not enter the objective
or optimizer. The unit test compares transformation Jacobians with ADMB
automatic derivatives in the regular and hybrid bounded modes.
