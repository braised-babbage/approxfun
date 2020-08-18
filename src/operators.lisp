(in-package :approxfun)

(defun lift-numeric-operator (op)
  "Lift an operator OP defined for numbers to one defined for functions."
  (lambda (&rest args)
    (let ((int nil))
      (dolist (f args)
	(when (chebyshev-approximant-p f)
	  (let ((int2 (chebyshev-approximant-interval f)))
	    (when (and int (not (interval= int int2)))
	      (error "Argument domain mismatch: ~A and ~A"
		     int int2))
	    (setf int int2))))
      (flet ((op-value (x)
               (apply op (mapcar (lambda (f) (function-value f x)) args))))
	(approxfun #'op-value :interval int)))))


(setf (symbol-function 'c+) (lift-numeric-operator #'+))
(setf (symbol-function 'c-) (lift-numeric-operator #'-))
(setf (symbol-function 'c*) (lift-numeric-operator #'*))
(setf (symbol-function 'c/) (lift-numeric-operator #'/))


(defun definite-integral (apfun)
  "Evaluate the definite integral of APFUN over its domain of definition."
  (let ((coeffs (chebyshev-approximant-coeffs apfun))
	(int (chebyshev-approximant-interval apfun)))
    ;; apfun(x)  = c(0) T_0(x) + c(1) T_1(x) + ...
    ;; \int_{-1}^1 T_k(x) = 0 if k odd, 2/(1-k^2) if k even
    (* (length-distortion -1d0 1d0
			  (interval-lower int) (interval-upper int))
       (realpart
	(loop :for i :from 0 :below (length coeffs)
              :for c := (aref coeffs i)
              :summing (if (oddp i)
                           0
                           (/ (* 2 c)
                              (- 1 (* i i)))))))))

(defun integrate (apfun &optional (constant 0))
  "Compute the indefinite integral of APFUN, where CONSTANT indicates the value
taken at the left endpoint of the domain."
  (let* ((coeffs (chebyshev-approximant-coeffs apfun))
	 (int (chebyshev-approximant-interval apfun))
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
    (let ((dx (length-distortion -1d0 1d0
				   (interval-lower int) (interval-upper int))))
	(unless (double= dx 1d0)
	  (loop :for i :from 0 :below (length coeffs)
		:do (setf (aref integrated i) (* (aref integrated i) dx)))))
    (let ((result (chebyshev-polynomial integrated :interval int)))
      ;; shift so that f(-1) = constant
      ;; TODO: this could be faster!!
      (c+ result
          (- constant (function-value result (interval-lower int)))))))

(defun differentiate (apfun)
  "Compute the derivative of APFUN."
  (let* ((int (chebyshev-approximant-interval apfun))
	 (coeffs (chebyshev-approximant-coeffs apfun))
         (n (length coeffs)))
    (if (= 1 n)
        (approxfun (constantly 0) :interval int)
	(let ((dcoeffs (make-array (1+ n) :element-type '(complex double-float))))
	  ;; FOR the interval [-1,1]
          ;; if p(x) = \sum_{k=0}^n a(k) T_k(x)
          ;; then p'(x) = \sum_{k=0}^{n-1} b(k) T_k(x)
          ;; with b(k-1) = b(k+1) + 2*k*a(k) for 2 <= k <= n
          ;; and b(n) = b(n+1) = 0, b(0) = b(2)/2 + a(1)
	  ;;
	  ;; In general, we really want f'(x) = g'(S(x))S'(x)
	  ;; where g is on [-1,1] and S maps int to [-1,1]
          (setf (aref dcoeffs n) #C(0d0 0d0)
                (aref dcoeffs (1- n)) #C(0d0 0d0))
          (loop :for i :from (1- n) :downto 2
                :do (setf (aref dcoeffs (1- i))
                          (+ (aref dcoeffs (1+ i))
			     (* 2d0 i (aref coeffs i)))))
          (setf (aref dcoeffs 0)
                (+ (/ (aref dcoeffs 2) 2d0)
                   (aref coeffs 1)))
	  (let ((dx (length-distortion (interval-lower int) (interval-upper int)
					 -1d0 1d0)))
	      (unless (double= dx 1d0)
		(loop :for i :from 0 :to n
		      :do (setf (aref dcoeffs i) (* (aref dcoeffs i) dx)))))
          (chebyshev-polynomial dcoeffs :interval int)))))

;;; TODO: mean, variance, std, norm
