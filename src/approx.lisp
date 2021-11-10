(in-package :approxfun)

(defstruct chebyshev-approximant
  "An approximation to a function via samples at Chebyshev points."
  name
  values
  coeffs
  interp-fn
  interval)

(defmethod domain ((ap chebyshev-approximant))
  (chebyshev-approximant-interval ap))

(defmethod print-object ((object chebyshev-approximant) stream)
  (print-unreadable-object (object stream :type t :identity t)
    (format stream "(~D)" (length (chebyshev-approximant-coeffs object)))
    (when (chebyshev-approximant-name object)
      (format stream " ~S " (chebyshev-approximant-name object)))
    (let ((d (chebyshev-approximant-interval object)))
      (format stream " on [~,3F, ~,3F]" (interval-lower d) (interval-upper d)))))

(defun chebyshev-polynomial (coeffs &key name (interval *default-interval*))
  "Construct a Chebyshev approximant directly from its coefficients."
  (let ((samples (samples-from-coefficients coeffs)))
    (make-chebyshev-approximant :name name
                                :values samples
                                :coeffs coeffs
                                :interp-fn (chebyshev-interpolate samples :interval interval)
                                :interval interval)))

(defun randomized-equality-check (obj1 obj2 &key (trials 10) (interval *default-interval*))
  "Check that function-like objects OBJ1 and OBJ2 agree on a random set of points."
  (let ((transform
          (affine-transformation
           0 1
           (interval-lower interval) (interval-upper interval))))
    (every (lambda (x)
             (double= (function-value obj1 x)
                      (function-value obj2 x)))
           (loop :for i :below trials
                 :collect (funcall transform (random 1d0))))))

(defparameter *log-max-chebyshev-samples* 15
  "The logarithm (base 2) of the maximum number of Chebyshev points to sample at.")

(defparameter *aliasing-check-num-samples* 2
  "The number of random samples to test as a check for aliasing.")

(defgeneric sample-at-chebyshev-points (obj n &key domain)
  (:documentation "Sample OBJ at N Chebyshev points on DOMAIN.")
  (:method ((fn function) n &key domain)
    (unless (typep domain 'interval)
      (domain-error "Invalid domain ~A" domain))
    (sample-fn-at-chebyshev-points fn n :interval domain))
  (:method ((fn chebyshev-approximant) n &key domain)
    (unless (domain-subset-p domain (domain fn))
      (domain-error "Domain mismatch: got ~A but expected ~A" domain (domain fn)))
    (sample-at-chebyshev-points (chebyshev-approximant-interp-fn fn)
                                n :domain domain)))

(defgeneric stopping-condition (fn ap &key domain)
  (:documentation "Return T if AP is a sufficiently close approximation to FN for adaptive sampling to terminate.")
  (:method ((fn function) (ap chebyshev-approximant) &key domain)
    (randomized-equality-check
     fn ap
     :trials *aliasing-check-num-samples*
     :interval domain))
  (:method ((fn chebyshev-approximant) (ap chebyshev-approximant) &key domain)
    (stopping-condition (chebyshev-approximant-interp-fn fn) ap :domain domain)))

(defun approxfun (fn &key num-samples name (interval *default-interval*))
  "Construct an approximation to the provided function.

If NUM-SAMPLES is provided, then the approximation is made from this fixed
number of function samples. Otherwise, adaptive sampling is used."
  (flet ((construct (vals)
           (make-chebyshev-approximant :name name
                                       :values vals
                                       :coeffs (chebyshev-coefficients vals)
                                       :interp-fn (chebyshev-interpolate vals :interval interval)
                                       :interval interval)))    
    (cond ((vectorp fn)
           (construct fn))
          (num-samples
           (construct (sample-at-chebyshev-points fn num-samples :domain interval)))
          (t
           ;; Adaptive search: we check on grids of size 2^d + 1, until we either
           ;; find an adequate approximation or we hit *MAX-CHEBYSHEV-SAMPLES*
           (loop :for d :from 4 :to *log-max-chebyshev-samples*
                 :for n := (1+ (expt 2 d))
                 :for vals := (sample-at-chebyshev-points fn n :domain interval)
                 :for coeffs := (chebyshev-coefficients vals)
                 :for cutoff := (coefficient-cutoff coeffs)
                 :when cutoff
                   :do (let ((ap (approxfun fn :num-samples (1+ cutoff) :name name :interval interval)))
                         (when (stopping-condition fn ap :domain interval)
                           (return ap)))
                 :finally (return (construct vals)))))))

(defun function-value (obj x)
  "Get the value of OBJ at point X."
  (cond ((numberp obj) obj)
        ((functionp obj) (funcall obj x))
        ((chebyshev-approximant-p obj)
         (funcall (chebyshev-approximant-interp-fn obj) x))
        (t
         (error "Unable to compute value of ~A at ~A." obj x))))
