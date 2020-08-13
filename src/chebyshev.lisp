(in-package :approxfun)

(defun chebyshev-points (n)
  "Construct an array of N Chebyshev points on the interval [-1,1]. "
  (let ((points (make-array n :element-type 'double-float)))
    (loop :for i :from 0 :below n
          :do (setf (aref points i) (cos (/ (* i pi) (1- n)))))
    points))

(defun sample-at-chebyshev-points (fn num-samples)
  "Sample a function FN at NUM-SAMPLES Chebyshev points."
  (let ((values (make-array num-samples :element-type 'double-float))
        (pts (chebyshev-points num-samples)))
    (loop :for i :from 0 :below num-samples
          :do (setf (aref values i)
                    (funcall fn (aref pts i))))
    values))

(defun chebyshev-coefficients (samples)
  "Get the coefficients of the Chebyshev interpolant of given SAMPLES."
  ;; we reflect to get one period [s(0), s(1), ..., s(n), s(n-1), ..., s(1)],
  ;; compute the FFT of this, and then truncate the results
  ;; NOTE: this could be more efficient using one of the DCT variants
  (let* ((n (length samples))
         (rlen (* 2 (1- n)))
         (reflected (make-array (list rlen)
                                :element-type '(complex double-float)))
         (results (make-array (list n)
                              :element-type '(complex double-float))))
    ;; construct reflected array
    (loop :for i :from 0 :below n
          :do (let ((val (complex (aref samples i))))
                (setf (aref reflected i) val)
                (when (< 0 i (1- n))
                  (setf (aref reflected (- rlen i)) val))))
    (let ((transformed (fft reflected)))
      ;; fft computes an unscaled transform, so we divide by n-1
      ;; note we are only keeping the first n coeffs
      (loop :for i :from 0 :below n
            :do (setf (aref results i)
                      (/ (aref transformed i) (1- n))))
      ;; we also need to scale appropriately at the boundary
      (setf (aref results 0)
            (/ (aref results 0) 2)
            (aref results (1- n))
            (/ (aref results (1- n)) 2))
      results)))

(defun chebyshev-interpolate (samples)
  "Given samples at the Chebyshev points, return a function computing the interpolant at an arbitrary point."
  ;; This is the so-called "Barycentric Interpolation", of Salzer
  ;; cf. https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
  (let* ((n (length samples))
         (xs (chebyshev-points n)))
    (lambda (x)
      (let ((num 0d0)
            (denom 0d0))
        (loop :for i :from 0 :below (length samples)
              :for w := 1 :then (- w)
              ;; TODO: should we be checking up to FP precision below?
              ;; I think the main thing is just to rule out actual divide by zero,
              ;; but its worth considering more carefully.
              :when (= x (aref xs i))
                :do (return (aref samples i))
              :do (let ((coeff
                          (/ (if (< 0 i (1- n)) w (/ w 2))
                             (- x (aref xs i)))))
                    (incf num (* (aref samples i) coeff))
                    (incf denom coeff))
              :finally (return (/ num denom)))))))
