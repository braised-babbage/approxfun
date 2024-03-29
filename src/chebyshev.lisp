(in-package :approxfun.coremath)

(defparameter *double-float-tolerance* 1d-15
  "Tolerance for double float calculations.")

(defun double= (x y)
  "Checks equality up to *DOUBLE-FLOAT-TOLERANCE*."
  (declare (type double-float x)
           (type double-float y))
  (< (abs (- x y))
     *double-float-tolerance*))

(defun double<= (x y)
  "Checks inequality up to *DOUBLE-FLOAT-TOLERANCE*"
  (declare (type double-float x)
           (type double-float y))
  (<= x (+ y *double-float-tolerance*)))

(defun clamp (x min max)
  (cond ((< x min) min)
        ((> x max) max)
        (t x)))

(defstruct (interval (:constructor %make-interval))
  "A repesentation of an interval [LOWER, UPPER]."
  lower
  upper)

(defun interval (lower upper)
  "Construct the interval [LOWER, UPPER]."
  (unless (< lower upper)
    (error "Unable to construct interval [~A, ~A]. " lower upper))
  (%make-interval :lower (coerce lower 'double-float)
                  :upper (coerce upper 'double-float)))

(defun interval= (int1 int2)
  "Are the two intervals equal?"
  (and (double= (interval-lower int1) (interval-lower int2))
       (double= (interval-upper int1) (interval-upper int2))))

(defun in-domain-p (x domain)
  "Is X in the DOMAIN?"
  (<= (interval-lower domain)
      x
      (interval-upper domain)))

(defun domain-subset-p (d1 d2)
  "Returns T if domain D1 is a subset of D2."
  (and (<= (interval-lower d2) (interval-lower d1))
       (<= (interval-upper d1) (interval-upper d2))))

(defparameter *default-interval* (interval -1 1)
  "The default interval.")

(defun affine-transformation (a0 b0 a1 b1)
  "Construct an affine transformation from [a0,b0] to [a1,b1].

Returns two values: the transformation itself, and its derivative."
  (when (= a0 b0)
    (error "AFFINE-TRANSFORMATION expects a nontrivial initial interval."))
  (let ((m (/ (- a1 b1)
              (- a0 b0)))
        (c (/ (- (* a0 b1) (* a1 b0))
              (- a0 b0))))
    (values (lambda (x)
              (+ (* m x) c))
            m)))

(defun length-distortion (a0 b0 a1 b1)
  "Get the factor by which length is stretched in the affine transformation from [a0,b0] to [a1,b1]."
  (/ (- a1 b1)
     (- a0 b0)))

(defun chebyshev-points (n &key (interval *default-interval*))
  "Construct an array of N Chebyshev points on the given INTERVAL. "
  (unless (> n 1)
    (error "Unable to construct Chebyshev points on grid of size ~D" n))
  (let ((points (make-array n :element-type 'double-float))
        (transform (affine-transformation -1d0 1d0
                                          (interval-lower interval) (interval-upper interval))))
    (loop :with m := (1- n)
          :for i :from m :downto 0
          :for k :from (- m) :by 2
          :do (setf (aref points i)
                    (funcall transform (sin (/ (* pi k) (* 2 m))))))
    points))

(defun sample-fn-at-chebyshev-points (fn num-samples &key (interval *default-interval*))
  "Sample a function FN at NUM-SAMPLES Chebyshev points."
  (let ((values (make-array num-samples :element-type 'double-float))
        (pts (chebyshev-points num-samples :interval interval)))
    (loop :for i :from 0 :below num-samples
          :do (setf (aref values i)
                    (funcall fn (aref pts i))))
    values))


(defun chebyshev-coefficients (samples)
  "Get the coefficients of the Chebyshev interpolant of given SAMPLES."
  ;; we reflect to get one period [s(0), s(1), ..., s(n-1), s(n-2), ..., s(1)],
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


(defun samples-from-coefficients (coeffs)
  "Construct samples at Chebyshev points from an array of Chebyshev coefficients."
  ;; coeffs is length n
  ;; apply inverse fourier transform
  ;; extract real part of first n values
  (let* ((n (length coeffs))
         (extended (make-array (* 2 (1- n)) :element-type '(complex double-float))))
    ;; zero pad to length 2(n-1)
    (loop :for i :from 0 :below (* 2 (1- n))
          :do (setf (aref extended i)
                    (if (< i n) (aref coeffs i) #C(0d0 0d0))))
    (let ((inverted (fft extended :direction :backward)))
      (let ((results (make-array n :element-type 'double-float)))
        (loop :for i :from 0 :below n
              :do (setf (aref results i) (realpart (aref inverted i))))
        results))))

(defun chebyshev-interpolate (samples &key (interval *default-interval*))
  "Given samples at the Chebyshev points, return a function computing the interpolant at an arbitrary point."
  ;; This is the so-called "Barycentric Interpolation", of Salzer
  ;; cf. https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
  (let* ((n (length samples))
         (xs (chebyshev-points n))
         (transform (affine-transformation (interval-lower interval) (interval-upper interval)
                                           -1d0 1d0)))
    (lambda (x)
      (let ((num 0d0)
            (denom 0d0)
            (x (funcall transform x)))
        (loop :for i :from 0 :below n
              :for w := 1 :then (- w)
              :for xi := (aref xs i)
              :for xdiff := (- x xi)
              :for fi := (aref samples i)
              ;; TODO: should we be checking up to FP precision below?
              ;; I think the main thing is just to rule out actual divide by zero,
              ;; but its worth considering more carefully.
              :when (zerop xdiff)
                :do (return fi)
              :do (let ((coeff (/ (if (< 0 i (1- n)) w (/ w 2))
                                  xdiff)))
                    (incf num (* fi coeff))
                    (incf denom coeff))
              :finally (return (/ num denom)))))))


(defun coefficient-cutoff (coeffs)
  "Find a cutoff point for the Chebyshev coefficients COEFFS.

This returns the last index of COEFFS which is deemed significant, or NIL if the
series if no such index is found. The heuristic used is described in Aurentz and
Trefthen, 'Chopping a Chebyshev Series'.

The tolerance *DOUBLE-FLOAT-TOLERANCE* is a relative tolerance, used to detect when the decay of
Chebyshev coefficients is deemed to be negligible."
  (declare (optimize (debug 3)))
  (let* ((n (length coeffs))
         (max-abs (loop :for i :from 0 :below n
                        :maximizing (abs (aref coeffs i))))
         (envelope (make-array n :element-type 'double-float)))
    (cond ((< n 17) nil)
          ((>= *double-float-tolerance* 1) 1)
          ((= 0d0 max-abs) 1)
          (t
           ;; Construct monotonic envelope
           (loop :with m := 0d0
                 :for i :from (1- n) :downto 0
                 :do (setf m (max m (abs (aref coeffs i)))
                           (aref envelope i) (/ m max-abs)))
           ;; scan for a plateau
           (multiple-value-bind (plateau-idx j2)
               (loop :for j1 :from 1 :below n
                     :for j2 := (round (+ (* 1.25 j1)
                                          5))
                     :when (<= n j2)
                       :do (return-from coefficient-cutoff nil)
                     :when (= 0d0 (aref envelope j1))
                       :do (return (values (1- j1) j2))
                     :do (let* ((e1 (aref envelope j1))
                                (e2 (aref envelope j2))
                                (r (* 3 (- 1 (/ (log e1) (log *double-float-tolerance*))))))
                           (when (< r (/ e2 e1))
                             (return (values (1- j1) j2))))
                     :finally (return-from coefficient-cutoff nil))
             ;; fix cutoff at a point where envelope + an affine function
             ;; is minimal
             (cond ((= 0d0 (aref envelope plateau-idx))
                    plateau-idx)
                   (t
                    (let ((j3 (loop :for i :from 0 :below n
                                    :until (< (aref envelope i)
                                              (expt *double-float-tolerance* (/ 7 6)))
                                    :finally (return i))))
                      (when (<= j3 j2)
                        (setf j2 j3
                              (aref envelope j2) (expt *double-float-tolerance* (/ 7 6))))
                      (loop :with min := 1d0
                            :with idx := 0
                            :for i :from 0 :to j2
                            :for cc := (+ (log (aref envelope i) 10)
                                          (/ (* i -1/3 (log *double-float-tolerance* 10))
                                             j2))
                            :when (< cc min)
                              :do (setf min cc
                                        idx i)
                            :finally (return (max (1- idx) 1)))))))))))

(defun chebyshev-differentiation-matrix (n &optional (half-width 1d0))
  "Compute the NxN Chebyshev differentiation matrix on an interval with the indicated HALF-WIDTH."
  (let* ((xs (chebyshev-points n))
         (mat (magicl:empty (list n n) :type 'double-float))
         (m (1- n)))
    (flet ((entry (i j)
             (let ((ci (if (< 0 i m) 1 2))
                   (cj (if (< 0 j m) 1 2))
                   (xi (aref xs i))
                   (xj (aref xs j)))
               (if (= i j)
                   (cond ((= i 0) (/ (1+ (* 2 m m)) 6 half-width))
                         ((= i m) (/ (1+ (* 2 m m)) -6 half-width))
                         (t (/ (- xi)
                               (* 2 (- 1 (* xj xj)))
                               half-width)))
                   (/ (* ci (if (evenp (+ i j)) 1 -1))
                      (* cj (- xi xj))
                      half-width)))))
      (magicl::into! #'entry mat))
    mat))
