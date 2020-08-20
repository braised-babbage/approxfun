(in-package #:approxfun)

(defun colleague-matrix (coeffs)
  (let ((n (loop :for i :from (1- (length coeffs)) :downto 0
		 :for c := (aref coeffs i)
		 :until (not (and (double= 0d0 (realpart c))
				  (double= 0d0 (imagpart c))))
		 :finally (return i))))
    (assert (> n 1))
    (let* ((mat (magicl:zeros (list n n) :type '(complex double-float))))
      ;; fill sub/superdiagonals
      (loop :for i :from 1 :below n
	    :do (setf (magicl:tref mat i     (1- i)) #C(0.5d0 0d0)
		      (magicl:tref mat (1- i)     i) #C(0.5d0 0d0)))
      (setf (magicl:tref mat 0 1) #C(1d0 0d0))
      ;; fill last row
      (let ((scale (/ 1 (* 2 (aref coeffs n)))))
	(loop :for j :from 0 :below n
	      :do (decf (magicl:tref mat (1- n) j)
			(* scale (aref coeffs j)))))
      mat)))


(defun roots (apfun)
  "Compute the roots of APFUN."
  (let ((int (chebyshev-approximant-interval apfun))
	(mat (colleague-matrix (chebyshev-approximant-coeffs apfun))))
    ;; roots are eigenvalues of the colleague matrix
    (let ((*double-float-tolerance* 1d-14))
      (loop :for c :in (magicl:eig mat)
	     :for re := (realpart c)
	     :for im := (imagpart c)
	     :when (and (double= 0d0 im)
			(double<= (interval-lower int) re)
			(double<= re (interval-upper int))) 
	       :collect (clamp re (interval-lower int) (interval-upper int))))))
