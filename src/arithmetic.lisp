(in-package #:approxfun)

(defmacro define-generic-binary-arithmetic (op impl-op &key documentation)
  (let ((binary-op-name (intern (format nil "TWO-ARG-~A" op))))
    `(progn
       (defgeneric ,binary-op-name (a b)
         (:method (a b)
           (error "Unable to perform ~A on ~A and ~A" ',op a b))
         (:method ((a number) (b number))
           (,impl-op a b))
         (:method ((a chebyshev-approximant) (b number))
           (approxfun (lambda (x)
                        (,impl-op (function-value a x) b))
                      :interval (chebyshev-approximant-interval a)))
         (:method ((a number) (b chebyshev-approximant))
           (approxfun (lambda (x)
                        (,impl-op a (function-value b x)))
                      :interval (chebyshev-approximant-interval b)))
         (:method ((a chebyshev-approximant) (b chebyshev-approximant))
           (let ((int-a (chebyshev-approximant-interval a))
                 (int-b (chebyshev-approximant-interval b)))
             (unless (interval= int-a int-b)
               (domain-error "Domain mismatch ~A ~A") int-a int-b)
             (approxfun (lambda (x)
                          (,impl-op (function-value a x) (function-value b x)))
                        :interval int-a))))
       (defun ,op (&rest args)
         ,@(if documentation
               (list documentation)
               nil)
         (reduce #',binary-op-name args)))))

(define-generic-binary-arithmetic + cl:+
  :documentation "Add two or more {numbers, approxfuns}.")

(define-generic-binary-arithmetic - cl:-
  :documentation "Subtract two or more {numbers, approfuns}.")

(define-generic-binary-arithmetic * cl:*
  :documentation "Multiply two or more {numbers, approfuns}.")

(define-generic-binary-arithmetic / cl:/
  :documentation "Divide two or more {numbers, approfuns}.")

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


