(defpackage #:approxfun
  (:use #:cl)
  (:export

   ;; approx.lisp
   #:approxfun
   #:*log-max-chebyshev-samples*

   ;; operators.lisp
   #:c+
   #:c-
   #:c*
   #:c/
   #:definite-integral

   ;; chebyshev.lisp
   #:chebyshev-points
   #:chebyshev-coefficients
   #:chebyshev-interpolate
   #:sample-at-chebyshev-points
   #:samples-from-coefficients

   ;; plot.lisp
   #:plot
   #:plot-coefficients
   ))
