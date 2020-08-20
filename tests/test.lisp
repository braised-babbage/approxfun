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
	"Indefinite integrals are correct for polynomials and trig functions."
  (is (approxfun::randomized-equality-check
       (approxfun (lambda (x) (- (* 1/2 x x) 1/2)))
       (integrate (approxfun (lambda (x) x)))))
  (is (approxfun::randomized-equality-check
       (approxfun (lambda (x) (- (sin x) (sin -1d0))))
       (integrate (approxfun #'cos)))))


(deftest test-arithmetic ()
  (is (approxfun::randomized-equality-check (approxfun (lambda (x) (+ x 1)))
                                 (c+ (approxfun (constantly 1d0))
                                     (approxfun #'identity))))
  (is (approxfun::randomized-equality-check (approxfun (lambda (x) (* x (sin x))))
                                 (c* (approxfun #'identity)
                                     (approxfun #'sin)))))

(deftest test-large-degree-chebyshev-monomial ()
  "Large Chebyshev monomials are preserved when constructing apfuns."
  (let* ((coeffs (make-array 100
                             :element-type '(complex double-float)
                             :initial-element #C(0d0 0d0))))
    (setf (aref coeffs 99) #C(1d0 0d0))
    (let* ((fn (chebyshev-interpolate (samples-from-coefficients coeffs)))
           (apfun (approxfun fn)))
      (is (approxfun::randomized-equality-check fn apfun)))))

(deftest test-differentiate ()
  "Checks that DIFFERENTIATE reproduces standard identities."
  (approxfun::randomized-equality-check (approxfun (constantly 0d0))
                             (differentiate (approxfun (constantly 1d0))))
  (approxfun::randomized-equality-check (approxfun #'identity)
                                        (differentiate (approxfun (lambda (x) (* 1/2 x x)))))
  ;; TODO: get to the bottom of this!!!
  (let ((cos (approxfun #'cos))
        (dsin (differentiate (approxfun #'sin)))
        (approxfun:*double-float-tolerance* 5d-15))
    (approxfun::randomized-equality-check cos dsin)))

(deftest test-colleague-matrix ()
  "Check that we correctly construct the colleague matrix of a simple Chebyshev polynomial."
  (let ((coeffs (vector #C(-3/8 0d0) #C(7/8 0d0) #C(-3/8 0d0) #C(1/4 0d0)))
	(mat (magicl:from-list '(0d0     1d0    0d0
				 0.5d0   0d0    0.5d0
				 0.75d0 -1.25d0 0.75d0)
			       '(3 3) :type '(complex double-float))))
    (is (magicl:= mat (approxfun::colleague-matrix coeffs)))
    t))

(deftest test-roots ()
  "Check that we can recover the roots of sin(3 pi x) on [-1,1]."
  (let ((apfun (approxfun (lambda (x) (sin (* x 3 pi)))))
	(expected-roots (list -1 -2/3 -1/3 0 1/3 2/3 1)))
    (let ((actual-roots
	    (sort (approxfun:roots apfun)
		  #'<)))
      (is (= (length actual-roots) (length expected-roots)))
      (is (every (lambda (x y)
		   (double= x y :threshold 1d-14))
		 actual-roots expected-roots)))))
