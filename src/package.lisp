;;; we stash some math stuff in approxfun.coremath, because in
;;; approxfun proper we're shadowing normal math ops

(defpackage #:approxfun.coremath
  (:use #:cl)
  (:export
   #:*double-float-tolerance*
   #:double=
   #:double<=
   #:clamp
   #:interval
   #:interval-lower
   #:interval-upper
   #:interval=
   #:in-domain-p
   #:*default-interval*
   #:affine-transformation
   #:length-distortion
   #:chebyshev-points
   #:sample-at-chebyshev-points
   #:chebyshev-coefficients
   #:samples-from-coefficients
   #:chebyshev-interpolate
   #:coefficient-cutoff
   #:chebyshev-differentiation-matrix))

(defpackage #:approxfun
  (:use #:cl #:approxfun.coremath)
  (:export
   ;; domain.lisp
   #:make-domain
   #:domain

   ;; approx.lisp
   #:approxfun
   #:function-value
   #:*log-max-chebyshev-samples*

   ;; operator.lisp
   #:diag
   #:const-operator
   #:I

   ;; chebyshev.lisp
   #:interval

   ;; arithmetic.lisp
   #:+
   #:-
   #:*
   #:/
   #:@

   ;; standard-functions.lisp
   #:sin
   #:cos
   #:tan
   #:exp
   #:log
   #:sqrt
   #:sinh
   #:cosh
   #:tanh
   #:expt

   ;; calculus.lisp
   #:definite-integral
   #:integrate
   #:differentiate
   #:maximum
   #:minimum

   ;; chebyshev.lisp
   #:chebyshev-points
   #:chebyshev-coefficients
   #:chebyshev-interpolate
   #:sample-at-chebyshev-points
   #:samples-from-coefficients
   #:*double-float-tolerance*

   ;; roots.lisp
   #:roots

   ;; solve.lisp
   #:dirichlet-boundary
   #:dirichlet-boundary-left
   #:dirichlet-boundary-right
   #:solve

   ;; plot.lisp
   #:plot
   #:plot-coefficients
   )
  (:shadow #:+ #:- #:* #:/
           #:sin #:cos #:tan
           #:exp #:log #:sqrt
           #:sinh #:cosh #:tanh
           #:expt))
