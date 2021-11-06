(in-package #:approxfun)

(defstruct operator
  domain
  forward-op
  matrix-constructor)

(defmethod domain ((op operator))
  (operator-domain op))

;;; TODO: define boundary conditions

;;; Constructors

(defun diag (apfun)
  (let ((domain (domain apfun))
        (interp-fn (chebyshev-approximant-interp-fn apfun)))
    (flet ((forward-op (f)
             (* apfun f))
           (matrix-ctr (num-samples)
             (magicl:from-diag
              (map 'list #'identity
                   (sample-at-chebyshev-points interp-fn
                                               num-samples
                                               :interval domain)))))      
      (make-operator :domain domain
                     :forward-op #'forward-op
                     :matrix-constructor #'matrix-ctr))))

(defun const-operator (domain const)
  (make-operator :domain domain
                 :forward-op (lambda (f) (* f const))
                 :matrix-constructor
                 (lambda (n)
                   (magicl:eye (list n n) :value (coerce const 'double-float)))))

(defun I (domain)
  (const-operator domain 1d0))
