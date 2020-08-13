(defpackage #:approxfun
  (:use #:cl)
  (:export

   ;; approx.lisp
   #:approxfun
   #:*log-max-chebyshev-samples*
   #:c+
   #:c-
   #:c*
   #:c/

   ;; chebyshev.lisp
   #:chebyshev-points
   #:chebyshev-coefficients
   #:chebyshev-interpolate))
