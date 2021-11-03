(fiasco:define-test-package #:approxfun-tests
  (:use #:cl)
  (:local-nicknames (:ap :approxfun))
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
  (let ((pts (ap:chebyshev-points 10)))
    (is (= 10 (length pts)))
    (is (= 1d0  (aref pts 0)))
    (is (= -1d0 (aref pts 9)))))

(deftest test-chebyshev-coefficients-simple-poly ()
  "Computed Chebyshev coefficients agree with expectations for simple polynomials."
  (flet ((f (x)
           (* 2 x x)))
    (is (double= #(1 0 0 0)
                 (ap:chebyshev-coefficients #(1d0 1d0 1d0 1d0))))
    (is (double= #(0 1 0 0)
                 (ap:chebyshev-coefficients
                  (ap:chebyshev-points 4))))
    (is (double= #(1 0 1 0)
                 (ap:chebyshev-coefficients
                  (map 'vector #'f
                       (ap:chebyshev-points 4)))))))


(deftest test-chebyshev-interpolation ()
  "Interpolation from samples at Chebyshev points agrees to FP precision with function values."
  (flet ((f (x)
           (* 2 x x)))
    (is (double= (map 'vector #'f (ap:chebyshev-points 4))
                 (map 'vector (ap:chebyshev-interpolate
                               (map 'vector #'f (ap:chebyshev-points 4)))
                      (ap:chebyshev-points 4))))
    (is (double= (map 'vector #'f
                      (ap:chebyshev-points 100))
                 (map 'vector (ap:chebyshev-interpolate
                               (map 'vector #'f (ap:chebyshev-points 4)))
                      (ap:chebyshev-points 100))))))

(deftest test-definite-integration ()
  "Definite integrals are correct for a few basic functions."
  (is (double= 0 (ap:definite-integral (ap:approxfun #'identity))))
  (is (double= 2/3 (ap:definite-integral (ap:approxfun (lambda (x) (* x x))))))
  (let ((c (ap:approxfun (lambda (x) (cos (* pi x)))))
        (s (ap:approxfun (lambda (x) (sin (* pi x))))))
    (is (double= 0 (ap:definite-integral c)))
    (is (double= 0 (ap:definite-integral s)))
    (is (double= 1 (ap:definite-integral (ap:* c c))))
    (is (double= 1 (ap:definite-integral (ap:* s s))))
    (is (double= 2 (ap:definite-integral (ap:+ (ap:* c c) (ap:* s s)))))))

(deftest test-coefficients-to-samples ()
  "We can construct Chebyshev samples from Chebyshev coefficients. "
  (let ((sin (ap:approxfun #'sin))
        (weird (ap:approxfun (lambda (x) (+ 1 (* 2 x x (exp x)))))))
    (is (double= (ap::chebyshev-approximant-values sin)
                 (ap::samples-from-coefficients
                  (ap::chebyshev-approximant-coeffs sin))))
    (is (double= (ap::chebyshev-approximant-values weird)
                 (ap::samples-from-coefficients
                  (ap::chebyshev-approximant-coeffs weird))))))


(deftest test-indefinite-integral ()
  "Indefinite integrals are correct for polynomials and trig functions."
  (is (ap::randomized-equality-check
       (ap:approxfun (lambda (x) (- (* 1/2 x x) 1/2)))
       (ap:integrate (ap:approxfun (lambda (x) x)))))
  (is (ap::randomized-equality-check
       (ap:approxfun (lambda (x) (- (sin x) (sin -1d0))))
       (ap:integrate (ap:approxfun #'cos)))))


(deftest test-arithmetic ()
  (is (ap::randomized-equality-check (ap:approxfun (lambda (x) (+ x 1)))
                                     (ap:+ (ap:approxfun (constantly 1d0))
                                           (ap:approxfun #'identity))))
  (is (ap::randomized-equality-check (ap:approxfun (lambda (x) (* x (sin x))))
                                     (ap:* (ap:approxfun #'identity)
                                           (ap:approxfun #'sin)))))

(deftest test-large-degree-chebyshev-monomial ()
  "Large Chebyshev monomials are preserved when constructing apfuns."
  (let* ((coeffs (make-array 100
                             :element-type '(complex double-float)
                             :initial-element #C(0d0 0d0))))
    (setf (aref coeffs 99) #C(1d0 0d0))
    (let* ((fn (ap:chebyshev-interpolate (ap:samples-from-coefficients coeffs)))
           (apfun (ap:approxfun fn)))
      (is (ap::randomized-equality-check fn apfun)))))

(deftest test-differentiate ()
  "Checks that DIFFERENTIATE reproduces standard identities."
  (ap::randomized-equality-check (ap:approxfun (constantly 0d0))
                                 (ap:differentiate (ap:approxfun (constantly 1d0))))
  (ap::randomized-equality-check (ap:approxfun #'identity)
                                 (ap:differentiate (ap:approxfun (lambda (x) (* 1/2 x x)))))
  ;; TODO: get to the bottom of this!!!
  (let ((cos (ap:approxfun #'cos))
        (dsin (ap:differentiate (ap:approxfun #'sin)))
        (ap:*double-float-tolerance* 5d-15))
    (is (ap::randomized-equality-check cos dsin))))

(deftest test-colleague-matrix ()
  "Check that we correctly construct the colleague matrix of a simple Chebyshev polynomial."
  (let ((coeffs (vector #C(-3/8 0d0) #C(7/8 0d0) #C(-3/8 0d0) #C(1/4 0d0)))
        (mat (magicl:from-list '(0d0     1d0    0d0
                                 0.5d0   0d0    0.5d0
                                 0.75d0 -1.25d0 0.75d0)
                               '(3 3) :type '(complex double-float))))
    (is (magicl:= mat (ap::colleague-matrix coeffs)))
    t))

(deftest test-roots ()
  "Check that we can recover the roots of sin(3 pi x) on [-1,1]."
  (let ((apfun (ap:approxfun (lambda (x) (sin (* x 3 pi)))))
        (expected-roots (list -1 -2/3 -1/3 0 1/3 2/3 1)))
    (let ((actual-roots
            (sort (ap:roots apfun)
                  #'<)))
      (is (= (length actual-roots) (length expected-roots)))
      (is (every (lambda (x y)
                   (double= x y :threshold 1d-14))
                 actual-roots expected-roots)))))

(deftest test-roots-insensitive-to-recursion-threshold ()
  "ROOTS gives the same results regardless of whether we trigger recursion."
  (let* ((apfun (ap:approxfun (lambda (x)
                                (* (sin (* 10 pi x))
                                   (exp (- (* x x)))))))
         (ncoeffs (length (ap::chebyshev-approximant-coeffs apfun))))
    (let ((default-roots (ap:roots apfun)))
      (is (= 21 (length default-roots)))
      (is (null
           (set-exclusive-or
            default-roots
            (ap:roots apfun :recursion-points-threshold (ceiling (/ ncoeffs 3)))
            :test (lambda (x y) (double= x y :threshold 1d-14))))))))

(deftest test-standard-functions ()
  "Irrational & transcendental functions apply to approxfuns correctly."
  (let ((x (ap:approxfun #'identity :interval (ap:interval 0.1 0.9)))
        (xs '(0.3d0 0.5d0 0.7d0)))
    (loop 
      :for fn :in '(sin cos tan exp log sqrt sinh cosh tanh)
      :for apfn := (funcall (intern (symbol-name fn) :ap)
                            x)
      :do (loop :for x0 :in xs
                :do (is (double= (funcall fn x0)
                                 (ap:@ apfn x0)))))
    (dolist (x0 xs)      
      (is (double= (expt x0 2)
                   (ap:@ (ap:expt x 2) x0))))))
