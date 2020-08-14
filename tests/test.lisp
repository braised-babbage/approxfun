(fiasco:define-test-package #:approxfun-tests
    (:use #:cl #:approxfun)
  (:export #:run-approxfun-tests))

(in-package #:approxfun-tests)

(defun run-approxfun-tests ()
  (run-package-tests :package ':approxfun-tests
                     :verbose t
                     :describe-failures t
                     :interactive nil))

(defun double= (obj1 obj2 &key (threshold 1d-15))
  (cond ((and (vectorp obj1) (vectorp obj2))
         (loop :for v1 :across obj1
               :for v2 :across obj2
               :when (> (abs (- v1 v2)) threshold)
                 :do (return-from double= nil)
               :finally (return t)))
        (t
         (< (abs (- obj1 obj2))
            threshold))))

(defun randomized-equality-check (obj1 obj2 &optional (trials 10))
  (every (lambda (x)
           (double= (approxfun::function-value obj1 x)
                    (approxfun::function-value obj2 x)))
         (loop :for i :below trials
               :collect (1- (random 1d0)))))

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

(deftest test-definite-integration ()
  "Definite integrals are correct for a few basic functions."
  (is (double= 0 (definite-integral (approxfun #'identity))))
  (is (double= 2/3 (definite-integral (approxfun (lambda (x) (* x x))))))
  (let ((c (approxfun (lambda (x) (cos (* pi x)))))
        (s (approxfun (lambda (x) (sin (* pi x))))))
    (is (double= 0 (definite-integral c)))
    (is (double= 0 (definite-integral s)))
    (is (double= 1 (definite-integral (c* c c))))
    (is (double= 1 (definite-integral (c* s s))))
    (is (double= 2 (definite-integral (c+ (c* c c) (c* s s)))))))

(deftest test-coefficients-to-samples ()
  "We can construct Chebyshev samples from Chebyshev coefficients. "
  (let ((sin (approxfun #'sin))
        (weird (approxfun (lambda (x) (+ 1 (* 2 x x (exp x)))))))
    (is (double= (approxfun::chebyshev-approximant-values sin)
                 (approxfun::samples-from-coefficients
                  (approxfun::chebyshev-approximant-coeffs sin))))
    (is (double= (approxfun::chebyshev-approximant-values weird)
                 (approxfun::samples-from-coefficients
                  (approxfun::chebyshev-approximant-coeffs weird))))))


(deftest test-indefinite-integral ()
  (is (randomized-equality-check
       (approxfun (lambda (x) (- (* 1/2 x x) 1/2)))
       (integrate (approxfun (lambda (x) x)))))
  (is (randomized-equality-check
       (approxfun (lambda (x) (- (sin x) (sin -1d0))))
       (integrate (approxfun #'cos)))))
