(fiasco:define-test-package #:approxfun-tests
    (:use #:cl #:approxfun)
  (:export #:run-approxfun-tests))

(in-package #:approxfun-tests)

(defun run-approxfun-tests ()
  (run-package-tests :package ':approxfun-tests
                     :verbose t
                     :describe-failures t
                     :interactive nil))

(defun double= (vec1 vec2 &key (threshold 1d-15))
  (assert (= (length vec1) (length vec2)))
  (loop :for v1 :across vec1
        :for v2 :across vec2
        :when (> (abs (- v1 v2)) threshold)
          :do (return-from double= nil))
  t)

(deftest test-chebyshev-points-size-and-sort ()
  "The array of Chebyshev points has the right size and order. "
  (let ((pts (chebyshev-points 10)))
    (is (= 10 (length pts)))
    (is (= 1d0  (aref pts 0)))
    (is (= -1d0 (aref pts 9)))))

(deftest test-chebyshev-coefficients-simple-poly ()
  "Computed Chebyshev coefficients agree with expectations for simple polynomials."
  (flet ((f (x)
           (* 2 x x)))
    (is (double= #(1 0 0 0)
                 (chebyshev-coefficients #(1d0 1d0 1d0 1d0))))
    (is (double= #(0 1 0 0)
                 (chebyshev-coefficients
                  (chebyshev-points 4))))
    (is (double= #(1 0 1 0)
                 (chebyshev-coefficients
                  (map 'vector #'f
                       (chebyshev-points 4)))))))

(deftest test-chebyshev-interpolation ()
  "Interpolation from samples at Chebyshev points agrees to FP precision with function values."
  (flet ((f (x)
           (* 2 x x)))
    (is (double= (map 'vector #'f (chebyshev-points 4))
                 (map 'vector (chebyshev-interpolate
                               (map 'vector #'f (chebyshev-points 4)))
                          (chebyshev-points 4))))
    (is (double= (map 'vector #'f
                      (chebyshev-points 100))
                 (map 'vector (chebyshev-interpolate
                               (map 'vector #'f (chebyshev-points 4)))
                      (chebyshev-points 100))))))
