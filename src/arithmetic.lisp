(in-package #:approxfun)

(defmacro define-generic-binary-arithmetic (op ptwise-impl mat-impl)
  (let ((binary-op-name (intern (format nil "TWO-ARG-~A" op))))
    `(progn
       (defgeneric ,binary-op-name (a b)
         (:method (a b)
           (error "Unable to perform ~A on ~A and ~A" ',op a b))
         (:method ((a number) (b number))
           (,ptwise-impl a b))
         (:method ((a chebyshev-approximant) (b number))
           (approxfun (lambda (x)
                        (,ptwise-impl (function-value a x) b))
                      :interval (chebyshev-approximant-interval a)))
         (:method ((a number) (b chebyshev-approximant))
           (approxfun (lambda (x)
                        (,ptwise-impl a (function-value b x)))
                      :interval (chebyshev-approximant-interval b)))
         (:method ((a chebyshev-approximant) (b chebyshev-approximant))
           (let ((int-a (chebyshev-approximant-interval a))
                 (int-b (chebyshev-approximant-interval b)))
             (unless (interval= int-a int-b)
               (domain-error "Domain mismatch ~A ~A") int-a int-b)
             (approxfun (lambda (x)
                          (,ptwise-impl (function-value a x) (function-value b x)))
                        :interval int-a)))
         ,@(when mat-impl
             `((:method ((a operator) (b number))
                 (,binary-op-name a (const-operator (domain a) b)))
               (:method ((a number) (b operator))
                 (,binary-op-name (const-operator (domain b) a) b))
               (:method ((a operator) (b chebyshev-approximant))
                 (,binary-op-name a (diag b)))
               (:method ((a chebyshev-approximant) (b operator))
                 (,binary-op-name (diag a) b))
               (:method ((a operator) (b operator))
                 (unless (domain= (domain a) (domain b))
                   (domain-error "Domain mismatch ~A ~A" (domain a) (domain b)))
                 (make-operator
                  :domain (domain a)
                  :forward-op (lambda (ap)
                                (,op (funcall (operator-forward-op a) ap)
                                     (funcall (operator-forward-op b) ap)))
                  :matrix-constructor (lambda (n)
                                        (,mat-impl (funcall (operator-matrix-constructor a) n)
                                                   (funcall (operator-matrix-constructor b) n)))))))))))

(define-generic-binary-arithmetic + cl:+ magicl:.+)

(defun + (&rest args)
  "Return the sum of its argments. With no args, returns 0."
  (if args
      (reduce #'two-arg-+ args)
      0))

(define-generic-binary-arithmetic - cl:- magicl:.-)

(defun - (a &rest args)
  "Subtract the second and all subsequent arguments from the first; or with one argument, negate the first argument."
  (if args
      (reduce #'two-arg-- (cons a args))
      (two-arg-- 0 a)))

(define-generic-binary-arithmetic * cl:* nil)

(defun * (&rest args)
  "Return the product of its arguments. With no args, returns 1."
  (if args
      (reduce #'two-arg-* args)
      1))

(defmethod two-arg-* ((a operator) (b number))
  (make-operator
   :domain (domain a)
   :forward-op (lambda (ap)
                 (* (funcall (operator-forward-op a) ap) b))
   :matrix-constructor (lambda (n)
                         (magicl:scale
                          (funcall (operator-matrix-constructor a) n)
                          b))))

(defmethod two-arg-* ((a number) (b operator))
  (two-arg-* b a))

(defmethod two-arg-* ((a operator) (b chebyshev-approximant))
  (two-arg-* a (diag b)))

(defmethod two-arg-* ((a chebyshev-approximant) (b operator))
  (two-arg-* b a))

(defmethod two-arg-* ((a operator) (b operator))
  (unless (domain= (domain a) (domain b))
    (domain-error "Domain mismatch ~A ~A" (domain a) (domain b)))
  (make-operator
   :domain (domain a)
   :forward-op (alexandria:compose (operator-forward-op a) (operator-forward-op b))
   :matrix-constructor (lambda (n)
                         (magicl:@ (funcall (operator-matrix-constructor a) n)
                                   (funcall (operator-matrix-constructor b) n)))))

(define-generic-binary-arithmetic / cl:/ nil)

(defun / (a &rest args)
  "Divide the first argument by each of the following arguments, in turn. With one argument, return reciprocal."
  (if args
      (reduce #'two-arg-/ (cons a args))
      (two-arg-/ 1 a)))

(defmethod two-arg-/ ((a operator) (b number))
  (two-arg-* a (cl:/ b)))

(defmethod two-arg-/ ((a number) (b operator))
  (error "Operator division is not defined. Consider a more explicit approach (e.g. using SOLVE)."))

(defmethod two-arg-/ ((a operator) (b chebyshev-approximant))
  (two-arg-* a (/ b)))

(defmethod two-arg-/ ((a chebyshev-approximant) (b operator))
  (error "Operator division is not defined. Consider a more explicit approach (e.g. using SOLVE)."))

(defmethod two-arg-/ ((a operator) (b operator))
  (error "Operator division is not defined. Consider a more explicit approach (e.g. using SOLVE)."))

(defgeneric @ (g f)
  (:documentation "Apply G to F.")
  (:method ((g chebyshev-approximant) (f real))
    (function-value g f))
  (:method ((g chebyshev-approximant) (f chebyshev-approximant))
    ;; this might fail due to bounds issues
    (approxfun (lambda (x) (function-value g (function-value f x)))
               :interval (chebyshev-approximant-interval f)))
  (:method ((op operator) (f chebyshev-approximant))
    (funcall (operator-forward-op op) f)))

(defgeneric solve (A b)
  (:documentation "Attempt to solve (@ A x) == b.")
  (:method ((A chebyshev-approximant) (b real))
    ;; we solve (- (@ A x) b) == 0
    (first (roots (- A b))))
  (:method ((A operator) (b chebyshev-approximant))
    ;; TODO: add adaptive soln
    (unless (domain= (domain A) (domain b))
      (domain-error "A has domain ~A, but b has domain ~A" (domain A) (domain b)))
    (let* ((n (length (chebyshev-approximant-values b)))
           (bvec (magicl::from-storage (chebyshev-approximant-values b) (list n 1)))
           (amat (funcall (operator-matrix-constructor A) n))
           (xvec (magicl:linear-solve amat bvec)))
      (approxfun (magicl::storage xvec) :interval (domain b)))))


