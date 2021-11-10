(in-package #:approxfun)

(defstruct operator
  domain
  forward-op
  matrix-constructor
  (derivative-order 0))

(defmethod domain ((op operator))
  (operator-domain op))

;;; Constructors

(defun diag (f)
  "Construct the operator g(x) |-> f(x)g(x) for a given approxfun F."
  (let ((domain (domain f))
        (interp-fn (chebyshev-approximant-interp-fn f)))
    (flet ((forward-op (g)
             (* f g))
           (matrix-ctr (n)
             (magicl:from-diag
              (map 'list #'identity
                   (sample-at-chebyshev-points interp-fn
                                               n
                                               :domain domain)))))
      (make-operator :domain domain
                     :forward-op #'forward-op
                     :matrix-constructor #'matrix-ctr))))

(defun const-operator (domain const)
  "Construct the operator g(x) |-> cg(x) for a given constant C."
  (make-operator :domain domain
                 :forward-op (lambda (f) (* f const))
                 :matrix-constructor
                 (lambda (n)
                   (magicl:eye (list n n) :value (coerce const 'double-float)))))

(defun I (domain)
  "Construct the identity operator on DOMAIN."
  (const-operator domain 1d0))

(defun D (domain)
  "Construct the derivative operator on DOMAIN."
  (make-operator :domain domain
                 :forward-op #'differentiate
                 :matrix-constructor #'chebyshev-differentiation-matrix
                 :derivative-order 1))

(defun D^2 (domain)
  "Construct the second-order derivative operator on DOMAIN."
  (* (D domain) (D domain)))

;;; construct matrix from operator

;; (defun operator-matrix (op domain n)
;;   (let ((elts (make-array n :element-type 'double-float :initial-element 0d0))
;;         (mat (magicl:zeros (list n n))))
;;     (dotimes (i n)
;;       (when (plusp i)
;;         (setf (aref elts (1- i)) 0d0))
;;       (setf (aref elts i) 1d0)
;;       (let* ((fn (approxfun elts :interval domain))
;;              (fn2 (funcall op fn))
;;              (res (chebyshev-approximant-values
;;                    (approxfun (chebyshev-approximant-interp-fn fn2) :num-samples n))))
;;         (dotimes (j n)
;;           (setf (magicl:tref mat j i) (aref res j)))))
;;     mat))
