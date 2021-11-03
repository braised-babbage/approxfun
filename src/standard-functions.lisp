(in-package #:approxfun)

(defmacro define-arithmetic-function (fn impl &key documentation domain)
  `(let ((domain ,domain))
     (declare (ignorable domain))
     (defun ,fn (x)
       ,@(if documentation (list documentation) nil)
       (typecase x
         (number (,impl x))
         (chebyshev-approximant
          (approxfun (lambda (y)
                       #|,@(when domain
                           (list `(unless (in-domain-p x domain)
                                    (domain-error "~A out of bounds for domain ~A" x domain))))|#
                       (funcall ',impl (function-value x y)))
                     :interval (chebyshev-approximant-interval x)))
         (otherwise
          (error 'type-error :datum x :expected-type `(or real chebyshev-approximant)))))))

(define-arithmetic-function sin cl:sin
  :documentation "Returns the sine of X.")

(define-arithmetic-function cos cl:cos
  :documentation "Returns the cosine of X.")

(define-arithmetic-function tan cl:tan
  ;; TODO: we might want to exclude pi/2, but for now it just gives us a huge number
  :documentation "Returns the tangent of X.")

(define-arithmetic-function exp cl:exp
  :documentation "Returns the exponential of X.")

(define-arithmetic-function log cl:log
  :domain (interval 0 most-positive-double-float) ; TODO we really want half-open interval
  :documentation "Returns the natural logarithm of X.")

(define-arithmetic-function sqrt cl:sqrt
  :domain (interval 0 most-positive-double-float)
  :documentation "Returns the non-negative square root of X.")

(define-arithmetic-function sinh cl:sinh
  :documentation "Returns the hyperbolic sine of X.")

(define-arithmetic-function cosh cl:cosh
  :documentation "Returns the hyperbolic cosine of X.")

(define-arithmetic-function tanh cl:tanh
  :documentation "Returns the hyperbolic tangent of X.")

(defun expt (base power)
  "Return BASE raised to the POWER."
  (etypecase base
    (number (cl:expt base power))
    (chebyshev-approximant
     (check-type power integer)
     (approxfun (lambda (x)
                  (cl:expt (function-value base x) power))
                :interval (chebyshev-approximant-interval base)))))
