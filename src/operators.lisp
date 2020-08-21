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
