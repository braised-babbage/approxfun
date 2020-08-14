(in-package :approxfun)

(defun lift-numeric-operator (op)
  "Lift an operator OP defined for numbers to one defined for functions."
  (lambda (&rest args)
    (flet ((op-value (x)
             (apply op (mapcar (lambda (f) (function-value f x)) args))))
      (approxfun #'op-value))))


(setf (symbol-function 'c+) (lift-numeric-operator #'+))
(setf (symbol-function 'c-) (lift-numeric-operator #'-))
(setf (symbol-function 'c*) (lift-numeric-operator #'*))
(setf (symbol-function 'c/) (lift-numeric-operator #'/))


(defun definite-integral (apfun)
  "Evaluate the definite integral of APFUN over its domain of definition."
  (let ((coeffs (chebyshev-approximant-coeffs apfun)))
    ;; apfun(x)  = c(0) T_0(x) + c(1) T_1(x) + ...
    ;; \int_{-1}^1 T_k(x) = 0 if k odd, 2/(1-k^2) if k even
    (realpart
     (loop :for i :from 0 :below (length coeffs)
           :for c := (aref coeffs i)
           :summing (if (oddp i)
                        0
                        (/ (* 2 c)
                           (- 1 (* i i))))))))

(defun integrate (apfun &optional (constant 0))
  "Compute the indefinite integral of APFUN, where CONSTANT indicates the value
taken at the left endpoint of the domain."
  (let* ((coeffs (chebyshev-approximant-coeffs apfun))
         (integrated (make-array (1+ (length coeffs))
                                 :element-type '(complex double-float)
                                 :initial-element #C(0d0 0d0))))
    ;; we have the identity
    ;; \int T_k(x) dx = T_{k+1}(x)/2(k+1) - T_{k-1}(x)/2(k-1)
    ;; with special cases for k = 0 and k = 1
    (setf (aref integrated 1) (aref coeffs 0))
    (when (< 1 (length coeffs))
      (setf (aref integrated 2) (/ (aref coeffs 1) 4)))
    (loop :for i :from 2 :below (length coeffs)
          :for c := (aref coeffs i)
          :do (incf (aref integrated (1+ i)) (/ c (* 2 (1+ i))))
          :do (decf (aref integrated (1- i)) (/ c (* 2 (1- i)))))
    (let ((result (chebyshev-polynomial integrated)))
      ;; shift so that f(-1) = constant
      (c+ result
          (- constant (function-value result -1d0))))))

;;; TODO: mean, variance, std, norm
