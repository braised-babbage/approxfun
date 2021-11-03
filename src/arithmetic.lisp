(in-package #:approxfun)

(define-condition domain-error (simple-error)
  ()
  (:documentation "An error associated with an invalid domain."))

(defun domain-error (format-control &rest format-args)
  "Signal a QUIL-EXPANSION-ERROR, incorporating information about the expansion context."
  (error 'domain-error
         :format-control format-control
         :format-arguments format-args))

(defmacro define-generic-binary-arithmetic (op &key documentation)
  (let ((binary-op-name (intern (format nil "TWO-ARG-~A" op)))
        (cl-op (intern (symbol-name op) :cl)))
    `(progn
       (defgeneric ,binary-op-name (a b)
         (:method (a b)
           (error "Unable to perform ~A on ~A and ~A" ',op a b))
         (:method ((a number) (b number))
           (,cl-op a b))
         (:method ((a chebyshev-approximant) (b number))
           (approxfun (lambda (x)
                        (,cl-op (function-value a x) b))
                      :interval (chebyshev-approximant-interval a)))
         (:method ((a number) (b chebyshev-approximant))
           (approxfun (lambda (x)
                        (,cl-op a (function-value b x)))
                      :interval (chebyshev-approximant-interval b)))
         (:method ((a chebyshev-approximant) (b chebyshev-approximant))
           (let ((int-a (chebyshev-approximant-interval a))
                 (int-b (chebyshev-approximant-interval b)))
             (unless (interval= int-a int-b)
               (domain-error "Domain mismatch ~A ~A") int-a int-b)
             (approxfun (lambda (x)
                         (,cl-op (function-value a x) (function-value b x)))
                       :interval int-a))))
       (defun ,op (&rest args)
         ,documentation
         (reduce #',binary-op-name args)))))

(define-generic-binary-arithmetic +
  :documentation "Add two or more {numbers, approfuns}.")

(define-generic-binary-arithmetic -
  :documentation "Add two or more {numbers, approfuns}.")

(define-generic-binary-arithmetic *
  :documentation "Add two or more {numbers, approfuns}.")

(define-generic-binary-arithmetic /
  :documentation "Add two or more {numbers, approfuns}.")

(defgeneric @ (g f)
  (:documentation "Apply G to F.")
  (:method ((g chebyshev-approximant) (f real))
    (function-value g f))
  (:method ((g chebyshev-approximant) (f chebyshev-approximant))
    ;; this might fail due to bounds issues
    (approxfun (lambda (x) (function-value g (function-value f x)))
               :interval (chebyshev-approximant-interval f))))

(defgeneric solve (A b)
  (:documentation "Attempt to solve (@ A x) == b.")
  (:method ((A chebyshev-approximant) (b real))
    ;; we solve (- (@ A x) b) == 0
    (first (roots (- A b)))))
