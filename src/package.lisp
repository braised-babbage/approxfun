(defpackage #:approxfun
  (:use #:cl)
  (:export

   ;; approx.lisp
   #:approxfun
   #:function-value
   #:*log-max-chebyshev-samples*

   ;; operators.lisp
   #:c+
   #:c-
   #:c*
   #:c/
   #:definite-integral
   #:integrate
   #:differentiate

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
   ))
