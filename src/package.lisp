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
   #:*default-interval*
   #:affine-transformation
   #:length-distortion
   #:chebyshev-points
   #:sample-at-chebyshev-points
   #:chebyshev-coefficients
   #:samples-from-coefficients
   #:chebyshev-interpolate
   #:coefficient-cutoff))

(defpackage #:approxfun
  (:use #:cl #:approxfun.coremath)
  (:export

   ;; approx.lisp
   #:approxfun
   #:function-value
   #:*log-max-chebyshev-samples*

   ;; arithmetic.lisp
   #:+
   #:-
   #:*
   #:/
   #:@
   #:solve

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

   ;; plot.lisp
   #:plot
   #:plot-coefficients
   )
  (:shadow #:+ #:- #:* #:/))
