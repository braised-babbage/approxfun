(in-package #:approxfun)

(defstruct (dirichlet-boundary
            (:constructor dirichlet-boundary (left right)))
  left
  right)

(define-condition incompatible-boundary-conditions (error)
  ((message :initarg :message
            :initform "Boundary conditions don't match operator order."
            :reader incompatible-boundary-conditions-message)
   (order :initarg :order
          :reader incompatible-boundary-conditions-order)
   (bc :initarg :bc
       :reader incompatible-boundary-conditions-bc))
  (:report (lambda (condition stream)
             (format stream "~A~&" (incompatible-boundary-conditions-message condition)))))

(defgeneric fixup-vector-for-boundary-conditions (vec bc order)
  (:documentation "Update (in-place) VEC to be compatible with the indicated boundary conditions.")
  (:method (vec (bc null) order)
    (declare (ignore bc order))
    vec)
  (:method (vec (bc dirichlet-boundary) order)
    (declare (ignore order))
    (let ((n (magicl:size vec)))
      (when (dirichlet-boundary-left bc)
        (setf (magicl:tref vec (1- n)) (dirichlet-boundary-left bc)))
      (when (dirichlet-boundary-right bc)
        (setf (magicl:tref vec 0) (dirichlet-boundary-right bc))))
    vec))

(defun fixup-matrix-dirichlet-left (mat)
  (let ((n (1- (magicl:nrows mat))))
    (setf (magicl:tref mat n n) 1d0)
    (loop :for j :below n
          :do (setf (magicl:tref mat n j) 0d0)))
  mat)

(defun fixup-matrix-dirichlet-right (mat)
  (setf (magicl:tref mat 0 0) 1d0)
  (loop :for j :from 1 :below (magicl:nrows mat)
        :do (setf (magicl:tref mat 0 j) 0d0))
  mat)

(defgeneric fixup-matrix-for-boundary-conditions (mat bc order)
  (:documentation "Update (in-place) MAT to be compatible with the indicated boundary conditions.")
  (:method (mat (bc null) order)
    (unless (zerop order)
      (error 'incompatible-boundary-conditions :order order :bc bc))
    mat)
  (:method (mat (bc dirichlet-boundary) order)
    (let ((left (dirichlet-boundary-left bc))
          (right (dirichlet-boundary-right bc)))
      (case order
        (0
         (when (or left right)
           (error 'incompatible-boundary-conditions :order order :bc bc))
         mat)
        (1
         (cond ((and left (null right))
                (fixup-matrix-dirichlet-left mat))
               ((and right (null left))
                (fixup-matrix-dirichlet-right mat))
               (t (error 'incompatible-boundary-conditions :order order :bc bc))))
        (2
         (unless (and left right)
           (error 'incompatible-boundary-conditions :order order :bc bc))
         (fixup-matrix-dirichlet-left mat)
         (fixup-matrix-dirichlet-right mat))
        (t
         (error "Solve does not support differential operators of order > 2, but got ~D." order))))))

(defgeneric solve (A b &key)
  (:documentation "Attempt to solve (@ A x) == b.")
  (:method ((A chebyshev-approximant) (b real) &key)
    ;; we solve (- (@ A x) b) == 0
    (first (roots (- A b))))
  (:method ((A operator) (b chebyshev-approximant) &key boundary-conditions)
    ;; TODO: add adaptive soln
    (unless (domain= (domain A) (domain b))
      (domain-error "A has domain ~A, but b has domain ~A" (domain A) (domain b)))
    (let ((order (operator-derivative-order A)))
      (when (> order 2)
        (error "Order > 2 differential operators currently unsupported in SOLVE."))
      (let* ((n (length (chebyshev-approximant-values b)))
             (bvec (magicl:from-array (chebyshev-approximant-values b) (list n)))
             (amat (funcall (operator-matrix-constructor A) n)))
        (fixup-vector-for-boundary-conditions bvec boundary-conditions order)
        (fixup-matrix-for-boundary-conditions amat boundary-conditions order)
        (let ((xvec (magicl:linear-solve amat bvec)))
          (approxfun (magicl::storage xvec) :interval (domain b)))))))