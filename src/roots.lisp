(in-package #:approxfun)

(defun colleague-matrix (coeffs)
  (let ((n (loop :for i :from (1- (length coeffs)) :downto 0
                 :for c := (aref coeffs i)
                 :until (not (and (double= 0d0 (realpart c))
                                  (double= 0d0 (imagpart c))))
                 :finally (return i))))
    (assert (> n 1))
    (let ((mat (magicl:zeros (list n n) :type '(complex double-float))))
      ;; fill sub/superdiagonals
      (loop :for i :from 1 :below n
            :do (setf (magicl:tref mat i     (1- i)) #C(0.5d0 0d0)
                      (magicl:tref mat (1- i)     i) #C(0.5d0 0d0)))
      (setf (magicl:tref mat 0 1) #C(1d0 0d0))
      ;; fill last row
      (let ((scale (cl:/ 1 (cl:* 2 (aref coeffs n)))))
        (loop :for j :from 0 :below n
              :do (decf (magicl:tref mat (1- n) j)
                        (cl:* scale (aref coeffs j)))))
      mat)))


(defun roots (apfun &key
                      (max-depth 10)
                      (recursion-points-threshold 100))
  "Compute the roots of APFUN."
  (labels ((find-roots (apfun depth)
             (let* ((lower (interval-lower (chebyshev-approximant-interval apfun)))
                    (upper (interval-upper (chebyshev-approximant-interval apfun)))
                    (mid (cl:/ (cl:+ lower upper) 2)))
               (cond ((< max-depth depth)
                      nil)
                     ((< recursion-points-threshold
                         (length (chebyshev-approximant-coeffs apfun)))
                      (union (find-roots (approxfun (chebyshev-approximant-interp-fn apfun)
                                                    :interval (interval lower mid))
                                         (1+ depth))
                             (find-roots (approxfun (chebyshev-approximant-interp-fn apfun)
                                                    :interval (interval mid upper))
                                         (1+ depth))
                             :test #'double=))
                     (t
                      (let ((mat (colleague-matrix (chebyshev-approximant-coeffs apfun)))
                            (transform (affine-transformation -1d0 1d0
                                                              lower upper)))
                        ;; roots are eigenvalues of the colleague matrix
                        (remove-duplicates
                         (loop :for c :in (magicl:eig mat)
                               :for x := (funcall transform (realpart c))
                               :for im := (imagpart c)
                               :when (and (double= 0d0 im)
                                          (double<= lower x)
                                          (double<= x upper))
                                 :collect x)
                         :test #'double=)))))))
    ;; fudge the fudge factor
    (let ((*double-float-tolerance* (cl:* 10 *double-float-tolerance*)))
      (find-roots apfun 0))))
