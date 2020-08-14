(in-package :approxfun)

(defstruct chebyshev-approximant
  "An approximation to a function via samples at Chebyshev points."
  name
  values
  coeffs
  interp-fn)

(defmethod print-object ((object chebyshev-approximant) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (when (chebyshev-approximant-name object)
      (format stream "~S" (chebyshev-approximant-name object)))))

(defun chebyshev-polynomial (coeffs &key name)
  "Construct a Chebyshev approximant directly from its coefficients."
  (let ((samples (samples-from-coefficients coeffs)))
    (make-chebyshev-approximant :name name
                                :values samples
                                :coeffs coeffs
                                :interp-fn (chebyshev-interpolate samples))))

(defparameter *log-max-chebyshev-samples* 15
  "The logarithm (base 2) of the maximum number of Chebyshev points to sample at.")

(defun approxfun (fn &key (num-samples nil) (name nil))
  "Construct an approximation to the provided function.

If NUM-SAMPLES is provided, then the approximation is made from this fixed
number of function samples. Otherwise, adaptive sampling is used."
  (if num-samples
        (let ((vals (sample-at-chebyshev-points fn num-samples)))
          (make-chebyshev-approximant :name name
                                      :values vals
                                      :coeffs (chebyshev-coefficients vals)
                                      :interp-fn (chebyshev-interpolate vals)))
        ;; Adaptive search: we check on grids of size 2^d + 1, until we either
        ;; find an adequate approximation or we hit *MAX-CHEBYSHEV-SAMPLES*
        (loop :for d :from 4 :to *log-max-chebyshev-samples*
              :for n := (1+ (expt 2 d))
              :for vals := (sample-at-chebyshev-points fn n)
              :for coeffs := (chebyshev-coefficients vals)
              :for cutoff := (coefficient-cutoff coeffs)
              :until cutoff
              :finally (return
                         (approxfun fn
                                    :num-samples (or (1+ cutoff) n)
                                    :name name)))))

(defun function-value (obj x)
  "Get the value of OBJ at point X."
  (cond ((numberp obj) obj)
        ((functionp obj) (funcall obj x))
        ((chebyshev-approximant-p obj)
         (funcall (chebyshev-approximant-interp-fn obj) x))
        (t
         (error "Unable to compute value of ~A at ~A." obj x))))
